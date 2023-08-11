/////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2023, The Regents of the University of California
// All rights reserved.
//
// BSD 3-Clause License
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

#pragma once

#include "db_sta/dbSta.hh"
#include "odb/db.h"
#include "odb/dbShape.h"
#include "utl/Logger.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <fstream>
#include <iostream>

namespace sec {

class secModule
{
 public:
  secModule();
  ~secModule();

  void init(odb::dbDatabase* db, sta::dbSta* sta, utl::Logger* logger);

  void evaluate();

  void readDesignMetrics(const std::string& filename);

  void readBaselineMetrics(const std::string& filename);

  void computeOverall();

 private:

  int64_t computeObstruction(uint width, uint length, odb::dbTechLayer* layer);

  // Evaluation Metrics
  // sec_ti_sts     = number of exploitable regions
  int64_t sec_ti_sts_ = 0;
  // sec_ti_sts_max = max number of sites across all exploitable regions
  int64_t sec_ti_sts_max_ = 0;
  // sec_ti_sts_med = median number of sites across all exploitable regions
  int64_t sec_ti_sts_med_ = 0;
  // sec_ti_sts_sum = total number of sites across all exploitable regions
  int64_t sec_ti_sts_sum_ = 0;
  // sec_ti_fts_sum = total number of free tracks across the whole layout and layers
  double sec_ti_fts_sum_ = 0;

  // Design Metrics -- Loaded from external scripts (file)
  bool design_loaded_ = false;
  double des_pwr_tot_ = 0;
  double des_ara_die_ = 0;
  double des_prf_wns_set_ = 0;
  double des_prf_wns_hld_ = 0;

  // Baseline Metrics
  bool baseline_loaded_ = false;
  int64_t sec_ti_sts_max_baseline_ = 0;
  int64_t sec_ti_sts_med_baseline_ = 0;
  int64_t sec_ti_sts_sum_baseline_ = 0;
  double sec_ti_fts_sum_baseline_ = 0;
  double des_pwr_tot_baseline_ = 0;
  double des_ara_die_baseline_ = 0;
  double des_prf_wns_set_baseline_ = 0;
  double des_prf_wns_hld_baseline_ = 0;

  // Threshold for a exploitable region
  int siteThreshold_ = 20; // Hardcoded 20 sites
  int siteConnection_ = 2; // Island must have at least 2 sites

  // Global
  sta::dbSta* sta_     = nullptr;
  odb::dbDatabase* db_ = nullptr;
  utl::Logger* logger_ = nullptr;
};

}  // namespace sec
