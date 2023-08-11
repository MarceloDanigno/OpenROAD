/////////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2023, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////

#include "sec/Sec.h"

#include "odb/db.h"
#include "utl/Logger.h"

namespace sec {

secModule::secModule() = default;

secModule::~secModule() = default;

void secModule::init(odb::dbDatabase* db, sta::dbSta* sta, utl::Logger* logger)
{
  db_ = db;
  logger_ = logger;
  sta_ = sta;

}

void secModule::evaluate()
{
  // TODO : do tests for empty chip, block, rows, insts

  odb::dbBlock* block = db_->getChip()->getBlock();
  odb::dbSet<odb::dbRow> rows = block->getRows();
  int siteWidth = block->getRows().begin()->getSite()->getWidth();
  int siteHeight = block->getRows().begin()->getSite()->getHeight();

  odb::dbSet<odb::dbNet> blockNets = block->getNets();
  odb::dbSet<odb::dbTrackGrid> trackGrids = block->getTrackGrids();

  // For track utilization -- Check FastRoute for other details
  int64_t routeTrackLength = 0;
  int64_t blockWirelength = 0;
  int64_t cellWirelength = 0;
  int blockWidth = block->getDieArea().dx();
  int blockHeight = block->getDieArea().dy();

  // Structures for finding exploitable regions
  std::map<int, std::vector<int>> rowPositions; // y [xmin, xmax]
  std::map<int, std::multiset<std::vector<int>>> cellPositions; // y [ {x1, x2}, {x1, x2} ... ]
  std::map<int, std::multiset<std::vector<int>>> islandPositions; // y [ {x1, x2}, {x1, x2} ... ]
  std::map<std::vector<int>, std::vector<std::vector<int>>> adjancencyMap; // {y, x1} [ {y, x1}, {y, x1} ... ]
  std::map<std::vector<int>, int> islandWidth; // {y, x1} width
  std::map<std::vector<int>, bool> visitedIslands; // {y, x1} 0/1
  // default for multisets = lesser sorting

  logger_->info(utl::SEC,
                1,
                "Running Security Evaluator on the design (site width = {} DBUs).",
                siteWidth);

  // Resetting previous evaluation
  sec_ti_sts_ = 0;
  sec_ti_sts_max_ = 0;
  sec_ti_sts_med_ = 0;
  sec_ti_sts_sum_ = 0;
  sec_ti_fts_sum_ = 0;
  std::multiset<int> exploitableRegionSet;

  // Number of DBUs that constitutes an exploitable region
  int siteConnection = siteWidth * siteConnection_; // hardcoded 2 * width DBU

  // Get the maximum track size for each layer
  // TODO : separate the fts/sts methods in functions
  for (odb::dbTrackGrid* grid : trackGrids){
    odb::dbTechLayer* layer = grid->getTechLayer();
    // Ignore routing on M1
    if(layer->getRoutingLevel() == 1){
      continue;
    }
    if (layer->getType() == odb::dbTechLayerType::ROUTING){
      int64_t trackLength = 0;
      if (layer->getDirection() == odb::dbTechLayerDir::HORIZONTAL){
        for (int i = 0; i < grid->getNumGridPatternsY(); i++){
          int offset, numTracks, spacing;
          grid->getGridPatternY(i, offset, numTracks, spacing);
          trackLength += (int64_t) blockWidth * numTracks;
        }
      } else {
        for (int i = 0; i < grid->getNumGridPatternsX(); i++){
          int offset, numTracks, spacing;
          grid->getGridPatternX(i, offset, numTracks, spacing);
          trackLength += (int64_t) blockHeight * numTracks;
        }
      }
      routeTrackLength += trackLength;
    }
  }

  logger_->info(utl::SEC,
              2,
              "Available Routing Track  = {} um.",
              ((double) routeTrackLength/4000));

  if (blockNets.empty()) {
    logger_->error(utl::SEC, 003, "Design with no nets.");
  }

  // Gets length of each metal segment for routing layers
  for (odb::dbNet* net : blockNets){
    uint wireCount = 0;
    uint viaCount = 0;
    net->getWireCount(wireCount, viaCount);
    if (wireCount == 0){
      continue;
    }
    // Supply/Special nets use SWires
    if (net->getSigType().isSupply() || net->isSpecial()){
      for (odb::dbSWire* swire : net->getSWires()){
        for (auto sbox : swire->getWires()) {
          // If the swire is a via, get the obstruction on routing layers
          if (sbox->isVia()){
            std::vector<odb::dbShape> viaBoxes;
            sbox->getViaBoxes(viaBoxes);
            for (odb::dbShape viaBox : viaBoxes){
              if (viaBox.getTechLayer()->getRoutingLevel() > 1) {
                if (viaBox.getTechLayer()->getDirection() == odb::dbTechLayerDir::HORIZONTAL){
                  blockWirelength += computeObstruction(viaBox.getDY(), viaBox.getDX(), viaBox.getTechLayer());
                } else {
                  blockWirelength += computeObstruction(viaBox.getDX(), viaBox.getDY(), viaBox.getTechLayer());
                }
              }
            }
          } else if (sbox->getTechLayer()->getRoutingLevel() > 1) {
            // Sbox already has a function for width and length
            blockWirelength += computeObstruction(sbox->getWidth(), sbox->getLength(), sbox->getTechLayer());
          }
        }
      }
    } else {
      // Signal nets use normal wires that can be divided into shapes
      odb::dbWire* wire = net->getWire();
      odb::dbWireShapeItr shapes;
      odb::dbShape s;
      for (shapes.begin(wire); shapes.next(s);) {
        // If the shape is a via, get the obstruction on routing layers
        if (s.isVia()){
          std::vector<odb::dbShape> viaBoxes;
          odb::dbShape::getViaBoxes(s, viaBoxes);
          for (odb::dbShape viaBox : viaBoxes){
            if (viaBox.getTechLayer()->getRoutingLevel() > 1) {
              if (viaBox.getTechLayer()->getDirection() == odb::dbTechLayerDir::HORIZONTAL){
                blockWirelength += computeObstruction(viaBox.getDY(), viaBox.getDX(), viaBox.getTechLayer());
              } else {
                blockWirelength += computeObstruction(viaBox.getDX(), viaBox.getDY(), viaBox.getTechLayer());
              }
            }
          }
        } else if (s.getTechLayer()->getRoutingLevel() > 1) {
          if (s.getTechLayer()->getDirection() == odb::dbTechLayerDir::HORIZONTAL){
            blockWirelength += computeObstruction(s.getDY(), s.getDX(), s.getTechLayer());
          } else {
            blockWirelength += computeObstruction(s.getDX(), s.getDY(), s.getTechLayer());
          }
        }
      }
    }
  }

  logger_->info(utl::SEC,
              4,
              "Nets and PDN Routing Utilization  = {}.",
              ((double) blockWirelength / routeTrackLength));

  sec_ti_fts_sum_ = ((double) blockWirelength / routeTrackLength) * 100;

  // Make a map for each row (y = key). Also initialize all the maps that depend on a row
  for (odb::dbRow* row : rows)
  {
    int xmin, xmax, y;
    odb::Rect bbox = row->getBBox();
    y = bbox.ul().getY();
    xmin = bbox.ul().getX();
    xmax = bbox.ur().getX();
    rowPositions[y] = {xmin, xmax};
    cellPositions[y] = {};
    islandPositions[y] = {};
  }

  // Make a map for the insts based on their y coordinate (y = key)
  // TODO : warn for inst out of row, inst no placed. Error if no instance is found on row
  for (odb::dbInst* inst : block->getInsts())
  {
    if (inst->getMaster()->isFiller() || inst->getMaster()->isEndCap()
     || inst->getMaster()->getType() == odb::dbMasterType::CORE_WELLTAP)
    {
      continue;
    }
    int xmin, xmax, y;
    odb::Rect bbox = inst->getBBox()->getBox();
    y = bbox.ul().getY();
    xmin = bbox.ul().getX();
    xmax = bbox.ur().getX();
    // TODO : check if y is in rowpositions before adding
    cellPositions[y].insert({xmin, xmax});
    odb::dbSet<odb::dbBox> obstructions = inst->getMaster()->getObstructions();
    for (odb::dbBox* obstruction : obstructions)
    {
      if (obstruction->getTechLayer()->getRoutingLevel() > 1) {
        if (obstruction->getTechLayer()->getDirection() == odb::dbTechLayerDir::HORIZONTAL){
          if (obstruction->getDY() > obstruction->getTechLayer()->getWidth()){
            cellWirelength += obstruction->getDX() + (obstruction->getDX() * (static_cast<int>(std::ceil(((double) obstruction->getDY()/2)/obstruction->getTechLayer()->getSpacing())) * 2));
          } else {
            cellWirelength += obstruction->getDX();
          }
        } else {
          if (obstruction->getDX() > obstruction->getTechLayer()->getWidth()){
            cellWirelength += obstruction->getDY() + (obstruction->getDY() * (static_cast<int>(std::ceil(((double) obstruction->getDX()/2)/obstruction->getTechLayer()->getSpacing())) * 2));
          } else {
            cellWirelength += obstruction->getDY();
          }
        }
      }
    }
  }

  logger_->info(utl::SEC,
              5,
              "Instances Routing Utilization  = {}.",
              ((double) cellWirelength / routeTrackLength));

  sec_ti_fts_sum_ = (1.0 - ((double) (cellWirelength + blockWirelength) / routeTrackLength)) * 100;

  // Populates islandPositions, a map of each island (2+ sites in the same row) based on their y coord
  for (auto row : rowPositions)
  {
    int y, xmin, xmax;
    y = row.first;
    xmin = row.second[0];
    xmax = row.second[1];
    if (cellPositions.find(y) != cellPositions.end())
    {
      int currentMin = xmin;

      // Check for exploitable regions between cells
      for (std::vector<int> cell : cellPositions[y])
      {
        int cellXmin, cellXmax;
        cellXmin = cell[0];
        cellXmax = cell[1];
        // Check if cell is out of bounds for the row
        if ((cellXmin < xmin) || (cellXmax > xmax))
        {
          continue;
        }
        else
        {
          int regionSize = cellXmin - currentMin;
          if (regionSize >= siteConnection)
          {
            islandPositions[y].insert({currentMin, cellXmin});
          }
          currentMin = cellXmax;
        }
      }

      // Compute between last cell and end of row
      int regionSize = xmax - currentMin;
      if (regionSize >= siteConnection)
      {
        islandPositions[y].insert({currentMin, xmax});
      }
    }
    else
    {
      // Row is empty
      int regionSize = xmax - xmin;
      if (regionSize >= siteConnection)
      {
        islandPositions[y].insert({xmin, xmax});
      }
    }
  }

  // Checks for adjacency between islands (only for the row below)
  for (auto& islandInterface : islandPositions)
  {
    int y;
    y = islandInterface.first;
    std::multiset<std::vector<int>> islandSet = islandInterface.second;
    for (auto& island : islandSet)
    {
      int xmin, xmax;
      xmin = island[0];
      xmax = island[1];

      if (adjancencyMap.find({y, xmin}) == adjancencyMap.end()){ // Add a blank connection
        adjancencyMap[{y, xmin}] = {};
      }

      islandWidth[{y, xmin}] = xmax - xmin;
      visitedIslands[{y, xmin}] = 0;

      if (islandPositions.find(y - siteHeight) != islandPositions.end())
      {
        std::multiset<std::vector<int>>& islandBelowSet = islandPositions[y - siteHeight];

        // TODO : reduce complexity by starting at the middle of the set, then checking left and right
        // TODO : keep doing this until you find all islands between xmin/xmax. Skip the rest

        // For each island, check if share a border that is >= than one site
        for (auto& belowIsland : islandBelowSet)
        {
          int belowXMin, belowXMax;
          bool connect = 0;
          belowXMin = belowIsland[0];
          belowXMax = belowIsland[1];

          if (belowXMin >= xmax) {
            break;
          }
          if (     (belowXMin <= xmin) && (belowXMax >= xmax))
          {
            connect = 1;
          }
          else if ((belowXMin <= xmin) && (belowXMax >= xmin) && ((belowXMax - xmin) >= (siteWidth)))
          {
            connect = 1;
          }
          else if ((belowXMin <= xmax) && (belowXMax >= xmax) && ((xmax - belowXMin) >= (siteWidth)))
          {
            connect = 1;
          }
          else if ((belowXMin >= xmin) && (belowXMax <= xmax))
          {
            connect = 1; // Since islands are already 2 sites wide, there is no need to check again
          }


          if (connect){

            adjancencyMap[{y, xmin}].push_back({y - siteHeight, belowXMin});
            if (adjancencyMap.find({y - siteHeight, belowXMin}) == adjancencyMap.end())
            {
              adjancencyMap[{y - siteHeight, belowXMin}] = {{y, xmin}};
            }
            else
            {
              adjancencyMap[{y - siteHeight, belowXMin}].push_back({y, xmin});
            }
          }
        }
      }
    }
  }

  // Check for regions that are higher than the threshold
  // Vector for debug print
  std::vector<std::vector<std::vector<int>>> exploitable_regions; // { {{a, b ,c} {a,b c} } }
  for (auto& islandConnection : adjancencyMap){
    if (visitedIslands[islandConnection.first]){
      continue;
    }

    std::vector<std::vector<int>> currentRegion;
    visitedIslands[islandConnection.first] = 1;
    std::vector<std::vector<int>> islandList = {islandConnection.second};
    int regionSiteSize = static_cast<int>(std::floor((double) islandWidth[islandConnection.first] / siteWidth));

    // For debug print
    int currY =  islandConnection.first[0];
    int currstartS =  islandConnection.first[1];
    int currendS =  islandWidth[islandConnection.first] + currstartS;
    std::vector<int> currentIsland = {currY, currstartS, currendS};
    currentRegion.push_back(currentIsland);

    while (!(islandList.empty())){
      std::vector<int> connectedIsland = islandList.back(); // Creates a copy of the last element
      islandList.pop_back();

      if (visitedIslands[connectedIsland]){
        continue;
      }

      visitedIslands[connectedIsland] = 1;

      // Keep adding adjacent island to the list until all of them are visited
      if (adjancencyMap.find(connectedIsland) != adjancencyMap.end()){
        for (auto& adjacentIsland : adjancencyMap[connectedIsland]){
          if (! visitedIslands[adjacentIsland]){
            islandList.push_back(adjacentIsland);
          }
        }
      }

      // Update the region size with the number of sites of each island
      regionSiteSize += static_cast<int>(std::floor((double) islandWidth[connectedIsland] / siteWidth));

      // For debug print
      int currY =  connectedIsland[0];
      int currstartS =  connectedIsland[1];
      int currendS =  islandWidth[connectedIsland] + currstartS;
      std::vector<int> currentConIsland = {currY, currstartS, currendS};
      currentRegion.push_back(currentConIsland);
    }

    if (regionSiteSize >= siteThreshold_)
    {
      // Is an exploitable region, update metrics
      exploitable_regions.push_back(currentRegion);
      sec_ti_sts_++;
      exploitableRegionSet.insert(regionSiteSize);
      sec_ti_sts_sum_ += regionSiteSize;
      if (regionSiteSize > sec_ti_sts_max_)
      {
        sec_ti_sts_max_ = regionSiteSize;
      }
    }
  }

  // Get median
  std::set<int>::iterator it = exploitableRegionSet.begin();
  std::advance(it, exploitableRegionSet.size() / 2);
  sec_ti_sts_med_ = *it;

  // Debug print
  if (logger_->debugCheck(utl::SEC, "secModule", 1)){
    for (int i = 0; i < exploitable_regions.size(); i++) {
          debugPrint(logger_,
                      utl::SEC,
                      "secModule",
                      1,
                      "  Region {} : {} islands.",
                      i,
                      exploitable_regions[i].size());
          for (auto& island : exploitable_regions[i]) {
              debugPrint(logger_,
                        utl::SEC,
                        "secModule",
                        1,
                        "      Island : {} row, {} start_site, {} end_site.",
                        island[0],
                        island[1],
                        island[2]);
          }
      }
  }

  const char *evaluationMessage =
    "Evaluation Results:\n\t\tsec_ti_sts = {}\n\t\tsec_ti_sts_max = {}"
    "\n\t\tsec_ti_sts_med = {}\n\t\tsec_ti_sts_sum = {}\n\t\tsec_ti_fts_sum = {}";

  logger_->info(utl::SEC,
            6,
            evaluationMessage,
            sec_ti_sts_, sec_ti_sts_max_, sec_ti_sts_med_, sec_ti_sts_sum_, sec_ti_fts_sum_);

}

int64_t secModule::computeObstruction(uint width, uint length, odb::dbTechLayer* layer){
  uint obstruction = length;
  if (width > layer->getWidth()){
    uint spacing = layer->getSpacing(width, length);
    if (spacing == 0){
      // M4 returns spacing == 0
      spacing = layer->getPitch() * 2;
    }
    obstruction += (length * (static_cast<int>(std::ceil((((double) width/2) + ((double) spacing))/layer->getPitch())) * 2));
  }
  return obstruction;
}

void secModule::readDesignMetrics(const std::string& filename){
  if (filename.empty()){
    logger_->error(utl::SEC, 7, "No design metrics file provided.");
  }
  // Initialize Metrics
  design_loaded_ = true;
  des_pwr_tot_ = 0;
  des_ara_die_ = 0;
  des_prf_wns_set_ = 0;
  des_prf_wns_hld_ = 0;

  // Read JSON input
  boost::property_tree::ptree jsonFile;
  std::ifstream ifs(filename, std::ios::binary);
  boost::property_tree::read_json(ifs, jsonFile);

  des_pwr_tot_ = jsonFile.get<double>("des_pwr_tot", 0);
  des_ara_die_ = jsonFile.get<double>("des_area_die", 0);
  des_prf_wns_set_ = jsonFile.get<double>("des_prf_WNS_set", 0);
  des_prf_wns_hld_ = jsonFile.get<double>("des_prf_WNS_hld", 0);
}

void secModule::readBaselineMetrics(const std::string& filename){
  if (filename.empty()){
    logger_->error(utl::SEC, 8, "No baseline metrics file provided.");
  }
  // Initialize Metrics
  baseline_loaded_ = true;
  sec_ti_sts_max_baseline_ = 0;
  sec_ti_sts_med_baseline_ = 0;
  sec_ti_sts_sum_baseline_ = 0;
  sec_ti_fts_sum_baseline_ = 0;
  des_pwr_tot_baseline_ = 0;
  des_ara_die_baseline_ = 0;
  des_prf_wns_set_baseline_ = 0;
  des_prf_wns_hld_baseline_ = 0;

  // Read JSON input
  boost::property_tree::ptree jsonFile;
  std::ifstream ifs(filename, std::ios::binary);
  boost::property_tree::read_json(ifs, jsonFile);

  sec_ti_sts_max_baseline_ = jsonFile.get<int>("sec_ti_sts_max", 0);
  sec_ti_sts_med_baseline_ = jsonFile.get<int>("sec_ti_sts_med", 0);
  sec_ti_sts_sum_baseline_ = jsonFile.get<int>("sec_ti_sts_sum", 0);
  sec_ti_fts_sum_baseline_ = jsonFile.get<double>("sec_ti_fts_sum", 0);
  des_pwr_tot_baseline_ = jsonFile.get<double>("des_pwr_tot", 0);
  des_ara_die_baseline_ = jsonFile.get<double>("des_area_die", 0);
  des_prf_wns_set_baseline_ = jsonFile.get<double>("des_prf_WNS_set", 0);
  des_prf_wns_hld_baseline_ = jsonFile.get<double>("des_prf_WNS_hld", 0);
}

void secModule::computeOverall(){
  if (!design_loaded_ || !baseline_loaded_){
    logger_->error(utl::SEC, 9, "Design or baseline metrics not loaded.");
  }

  double sec_ti_sts = ((1.0/2) * ((double) sec_ti_sts_sum_ / sec_ti_sts_sum_baseline_)) +
                      ((1.0/3) * ((double) sec_ti_sts_max_ / sec_ti_sts_max_baseline_)) +
                      ((1.0/6) * ((double) sec_ti_sts_med_ / sec_ti_sts_med_baseline_));
  double sec = ((1.0/2) * sec_ti_sts) + ((1.0/2) * (sec_ti_fts_sum_ / sec_ti_fts_sum_baseline_));
  double des_prf = ((1.0/2) * (des_prf_wns_set_ / des_prf_wns_set_baseline_)) +
                  ((1.0/2) * (des_prf_wns_hld_ / des_prf_wns_hld_baseline_));
  double des = ((1.0/3) * (des_pwr_tot_ / des_pwr_tot_baseline_)) +
               ((1.0/3) * (des_ara_die_ / des_ara_die_baseline_)) +
               ((1.0/3) * des_prf);
  double overall = ((1.0/2) * sec) + ((1.0/2) * des);

  logger_->info(utl::SEC,
            10,
            "Overall Security Score = {}.",
            overall);
}

}  // namespace sec
