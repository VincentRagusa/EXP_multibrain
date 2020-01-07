//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "GenomicMultiBrainWorld.h"
#include "../../Genome/CircularGenome/CircularGenome.h"

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::groupNamePL = Parameters::register_parameter("WORLD_TEST_NAMES-groupNameSpace", (std::string) "root::", "namespace of group to be evaluated");

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::genomeNamePL = Parameters::register_parameter( "WORLD_TEST_NAMES-genomeNameSpace", (std::string) "root::", "namespace for parameters used to define brain");

GenomicMultiBrainWorld::GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_):
  AbstractWorld(PT_) {
  // columns to be added to ave file
  popFileColumns.clear();
  popFileColumns.push_back("score");
}

void
GenomicMultiBrainWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug) {
  double score = 0.0;

  auto genome = std::dynamic_pointer_cast<CircularGenome<int>>(org->genomes[genomeNamePL->get(PT)]);
  if (genome == nullptr){
    std::cout << "GARBO set your genome site type to int ya fool (not really i'm the fool)" << std::endl;
    exit(1);
  }

  for (auto site:genome->sites){
    score += 100 - std::abs(100 - site);
  }
  
  org->dataMap.append("score", score);
  
}

void
GenomicMultiBrainWorld::evaluate(std::map<std::string, std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug) {
  for (auto& org:groups[groupNamePL->get(PT)]->population){
    evaluateSolo(org, analyze, visualize, debug);
  }
}

std::unordered_map<std::string, std::unordered_set<std::string>>
GenomicMultiBrainWorld::requiredGroups() {
  // return {{groupNamePL->get(PT), {"B:" + genomeNamePL->get(PT) + ",1," + std::to_string(numberOfOutputsPL->get(PT))}}};
  return {{groupNamePL->get(PT), {"G:" + genomeNamePL->get(PT)}}};
  // requires a root group and a brain (in root namespace) and no additional genome,
  // the brain must have 1 input, and the variable numberOfOutputs outputs
}



