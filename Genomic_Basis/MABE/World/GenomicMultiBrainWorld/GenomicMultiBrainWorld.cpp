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

std::shared_ptr<ParameterLink<int>>
GenomicMultiBrainWorld::modePL = Parameters::register_parameter("WORLD_TEST-mode", 0, "0 = bit outputs before adding, 1 = add outputs");

std::shared_ptr<ParameterLink<int>>
GenomicMultiBrainWorld::numberOfOutputsPL = Parameters::register_parameter("WORLD_TEST-numberOfOutputs", 10, "number of outputs in this world");

std::shared_ptr<ParameterLink<int>>
GenomicMultiBrainWorld::evaluationsPerGenerationPL = Parameters::register_parameter("WORLD_TEST-evaluationsPerGeneration", 1, "Number of times to test each Genome per generation (useful with non-deterministic brains)");

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::groupNamePL = Parameters::register_parameter("WORLD_TEST_NAMES-groupNameSpace", (std::string) "root::", "namespace of group to be evaluated");

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::brainNamePL = Parameters::register_parameter( "WORLD_TEST_NAMES-brainNameSpace", (std::string) "root::", "namespace for parameters used to define brain");

GenomicMultiBrainWorld::GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_):
  AbstractWorld(PT_) {
  // columns to be added to ave file
  popFileColumns.clear();
  popFileColumns.push_back("score");
  popFileColumns.push_back("score_VAR"); // specifies to also record the variance (performed automatically because _VAR)
}

void
GenomicMultiBrainWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug) {
  for (int r = 0; r < evaluationsPerGenerationPL->get(PT); r++) {
    double score = 0.0;
    org->dataMap.append("score", score);
    if (visualize) std::cout << "organism with ID " << org->ID << " scored " << score << std::endl;
  }
}

void
GenomicMultiBrainWorld::evaluate(std::map<std::string, std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug) {
  int popSize = groups[groupNamePL->get(PT)]->population.size();
  for (int i = 0; i < popSize; i++) {
    evaluateSolo(groups[groupNamePL->get(PT)]->population[i], analyze, visualize, debug);
  }
}

std::unordered_map<std::string, std::unordered_set<std::string>>
GenomicMultiBrainWorld::requiredGroups() {
  return {{groupNamePL->get(PT), {"B:" + brainNamePL->get(PT) + ",1," + std::to_string(numberOfOutputsPL->get(PT))}}};
  // requires a root group and a brain (in root namespace) and no additional genome,
  // the brain must have 1 input, and the variable numberOfOutputs outputs
}



