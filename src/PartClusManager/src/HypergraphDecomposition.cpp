////////////////////////////////////////////////////////////////////////////////
//// Authors: Mateus Fogaça, Isadora Oliveira and Marcelo Danigno
////
////          (Advisor: Ricardo Reis and Paulo Butzen)
////
//// BSD 3-Clause License
////
//// Copyright (c) 2020, Federal University of Rio Grande do Sul (UFRGS)
//// All rights reserved.
////
//// Redistribution and use in source and binary forms, with or without
//// modification, are permitted provided that the following conditions are met:
////
//// * Redistributions of source code must retain the above copyright notice, this
////   list of conditions and the following disclaimer.
////
//// * Redistributions in binary form must reproduce the above copyright notice,
////   this list of conditions and the following disclaimer in the documentation
////   and/or other materials provided with the distribution.
////
//// * Neither the name of the copyright holder nor the names of its
////   contributors may be used to endorse or promote products derived from
////   this software without specific prior written permission.
////
//// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//// POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////

#include "HypergraphDecomposition.h"
#include "opendb/db.h"
#include <fstream>

namespace PartClusManager{

void HypergraphDecomposition::init(int dbId){
	_db = odb::dbDatabase::getDatabase(dbId);
	_chip = _db->getChip();
	_block = _chip->getBlock();
}

void HypergraphDecomposition::constructMap(Hypergraph & hypergraph, unsigned maxVertexWeight){
	for (odb::dbNet* net : _block->getNets()){
		int nITerms = (net->getITerms()).size();
		int nBTerms = (net->getBTerms()).size();
		if (nITerms + nBTerms < 2)
			continue;
		for (odb::dbBTerm* bterm : net->getBTerms()){
			for (odb::dbBPin * pin : bterm->getBPins()){
				if (!hypergraph.isInMap(bterm->getName())){ 
					odb::dbBox * box = pin->getBox();
					long long int length = box->getLength();
					long long int width = box->getWidth();
					long long int area = length * width;
					int nextIdx = hypergraph.computeNextVertexIdx(); 
					hypergraph.addMapping(bterm->getName(), nextIdx);
					hypergraph.addVertexWeight(area);
				}
			}
		}
	
		for (odb::dbITerm * iterm : net->getITerms()){
			odb::dbInst * inst = iterm->getInst();
			if (!hypergraph.isInMap(inst->getName())){
				odb::dbBox * bbox = inst->getBBox();
				long long int length = bbox->getLength();
				long long int width = bbox->getWidth();
				long long int area = length * width;	
				int nextIdx = hypergraph.computeNextVertexIdx();
				hypergraph.addMapping(inst->getName(), nextIdx);
				hypergraph.addVertexWeight(area);
			}
		}
	}
	hypergraph.computeVertexWeightRange(maxVertexWeight);
}


void HypergraphDecomposition::createHypergraph(Hypergraph &hypergraph, std::vector<short> clusters , short currentCluster){
	for (odb::dbNet* net : _block->getNets()){
		int nITerms = (net->getITerms()).size();
		int nBTerms = (net->getBTerms()).size();
		if (nITerms + nBTerms < 2)
			continue;
		
		int driveIdx = -1;
		std::vector<int> netVertices;
		int nextPtr =  hypergraph.computeNextRowPtr();
		hypergraph.addRowPtr(nextPtr);
		hypergraph.addEdgeWeightNormalized(1);
		
		for (odb::dbBTerm* bterm : net->getBTerms()){
			for (odb::dbBPin * pin : bterm->getBPins()){
				int mapping = hypergraph.getMapping(bterm->getName());
				if (clusters[mapping] == currentCluster){
		                        if (bterm->getIoType() == odb::dbIoType::INPUT)
                		                driveIdx = mapping;
					else
						netVertices.push_back(mapping);

				} 
			}
		}
	
		for (odb::dbITerm * iterm : net->getITerms()){
			odb::dbInst * inst = iterm->getInst();
			int mapping = hypergraph.getMapping(inst->getName());
			if (clusters[mapping] == currentCluster){
				if (driveIdx == -1)
                        		if (iterm->isOutputSignal()){
                                		driveIdx = mapping;
						continue;
					}
				netVertices.push_back(mapping);
				
			} 
		}
		if (driveIdx != -1)
			hypergraph.addColIdx(driveIdx);
		for (int vertex : netVertices){
			hypergraph.addColIdx(vertex);}
	}
	int nextPtr =  hypergraph.computeNextRowPtr();
	hypergraph.addRowPtr(nextPtr);
}

GraphType HypergraphDecomposition::resolveModel(std::string graphModel){
	if (graphModel == "clique"){
		return CLIQUE;
	}
	if (graphModel == "star"){
		return STAR;
	}
	if (graphModel == "hybrid"){
		return HYBRID;
	}
}

void HypergraphDecomposition::toGraph(Hypergraph &hypergraph, Graph & graph, std::string graphModelS, unsigned weightingOption, 
				unsigned maxEdgeWeight, unsigned threshold){
	_weightingOption = weightingOption;
	GraphType graphModel = resolveModel(graphModelS);
	std::vector<int> colIdx = hypergraph.getColIdx();
	adjMatrix.resize(hypergraph.getNumVertex());
	for(int i=0; i < hypergraph.getNumRowPtr() -1; i++){
		std::vector<int> net;
		std::vector<int>::iterator begin = colIdx.begin() + hypergraph.getRowPtr(i);
		std::vector<int>::iterator end = colIdx.end();
		end = colIdx.begin() + hypergraph.getRowPtr(i+1);	
		net.assign(begin, end);

		switch(graphModel){
			case CLIQUE:
				if (net.size() > threshold)
					continue;
				createCliqueGraph(graph, net);
				break;
			case STAR:
				createStarGraph(graph, net);
				break;
			case HYBRID:
				if (net.size() > threshold){ 
					createStarGraph(graph, net);
				} else{ 
					createCliqueGraph(graph, net);
				}
				break;
		}
	
	}
	createCompressedMatrix(graph);
	graph.assignVertexWeight(hypergraph.getVertexWeight());
	graph.computeEdgeWeightRange(maxEdgeWeight);
}

float HypergraphDecomposition::computeWeight(int nPins){
	switch(_weightingOption){
		case 1:
			return 1.0/(nPins-1);
		case 2:
			return 4.0/(nPins*(nPins-1));
		case 3:
			return 4.0/(nPins*nPins - (nPins % 2));
		case 4:
			return 6.0/(nPins*(nPins+1));
		case 5:
			return pow((2.0/nPins),1.5);
		case 6:
			return pow((2.0/nPins),3);
		case 7:
			return 2.0/nPins;
	}
}

void HypergraphDecomposition::connectStarPins(int firstPin, int secondPin, float weight){
	bool isConnected = false;
	if (firstPin == secondPin)
		return;
        if (adjMatrix[firstPin].find(secondPin) != adjMatrix[firstPin].end())
			isConnected = true;
	if (!isConnected)
		adjMatrix[firstPin][secondPin] = weight;
	
}

void HypergraphDecomposition::connectPins(int firstPin, int secondPin, float weight){
	bool isConnected = false;
	if (firstPin == secondPin)
		return;
        if (adjMatrix[firstPin].find(secondPin) != adjMatrix[firstPin].end())
			isConnected = true;

	if (!isConnected){
		adjMatrix[firstPin][secondPin] = weight;
	} else{
		adjMatrix[firstPin][secondPin] += weight;
	}
	
}

void HypergraphDecomposition::createStarGraph(Graph & graph, std::vector<int> net){
	float weight = 1; 
	int driveIdx = 0;
	for (int i =1; i < net.size(); i++){
		connectStarPins(net[i], net[driveIdx], weight);
		connectStarPins(net[driveIdx], net[i], weight);
	}
}

void HypergraphDecomposition::createCliqueGraph(Graph & graph, std::vector<int> net){
	float weight = computeWeight(net.size()); 
	for (int i =0; i < net.size(); i++){
		for (int j = i+1; j < net.size(); j++){
			connectPins(net[i], net[j], weight);		
			connectPins(net[j], net[i], weight);		
		}
	}

}

void HypergraphDecomposition::createCompressedMatrix(Graph & graph){

	for (std::map<int,float> & node : adjMatrix){
		int nextPtr = graph.computeNextRowPtr();
		graph.addRowPtr(nextPtr);
		for (std::pair<const int,float> & connection : node){
			graph.addEdgeWeight(connection.second);
			graph.addColIdx(connection.first);
		}
	}	
	
}



}
