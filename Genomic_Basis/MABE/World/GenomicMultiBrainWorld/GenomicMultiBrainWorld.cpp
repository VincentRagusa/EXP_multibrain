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
GenomicMultiBrainWorld::groupNamePL = Parameters::register_parameter("WORLD_GMB_NAMES-groupNameSpace", (std::string) "root::", "namespace of group to be evaluated");

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::genomeNamePL = Parameters::register_parameter( "WORLD_GMB_NAMES-genomeNameSpace", (std::string) "root::", "namespace for parameters used to define brain");

std::shared_ptr<ParameterLink<int>>
GenomicMultiBrainWorld::gameCodePL = Parameters::register_parameter( "WORLD_GMB-gameCode", (int) 0, "it picks the game yo");


GenomicMultiBrainWorld::GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_):
  AbstractWorld(PT_) {
  gameCode = gameCodePL->get(PT);
  convertCSVListToVector(genomeNamePL->get(PT), genomeNames);

  // columns to be added to ave file
  popFileColumns.clear();
  popFileColumns.push_back("score");
  for (auto name : genomeNames) { 
		popFileColumns.push_back(name + "score");
	}
  if (gameCode == 2){
    popFileColumns.push_back("mean");
    popFileColumns.push_back("variance");
    for (auto name : genomeNames) { 
      popFileColumns.push_back(name + "mean");
      popFileColumns.push_back(name + "variance");
    }
  
  }
}

void
GenomicMultiBrainWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug) {
  for (auto name : genomeNames) {
    // if (Global::update == 0){
    //   //randomize genomes
    //   org->genomes[name]->fillRandom();  
    // }

    double score = 0.0;

    auto genome = std::dynamic_pointer_cast<CircularGenome<int>>(org->genomes[name]);
    if (genome == nullptr){
      std::cout << "GARBO set your genome site type to int ya fool (not really i'm the fool)" << std::endl;
      exit(1);
    }

    if (gameCode == 0){
      for (auto site:genome->sites){
        score += 100 - std::abs(100 - site); //100 is the target site value, assumed to be within the alphabet size
      }
    }
    else if (gameCode == 1){
      std::fill(site_symbol_counts.begin(), site_symbol_counts.end(), 0);
      for (auto site:genome->sites){ //TODO BUG fix the last site in the genome bug for this game
        ++site_symbol_counts[site];
      }
      auto max = *std::max_element(site_symbol_counts.begin(), site_symbol_counts.end());
      score = 5000 - max; //5000 is the assumed genome length
    }
    else if (gameCode == 2){
      auto mean = std::accumulate(genome->sites.begin(), genome->sites.end()-1, 0.0) / (genome->sites.size()-1); //last site is borked
      double variance = 0;
      variance = *std::max_element(genome->sites.begin(), genome->sites.end()-1) - *std::min_element(genome->sites.begin(), genome->sites.end()-1); //NOT REWAL
      // std::cout << variance << std::endl;
      // // for (auto site:genome->sites){
      // for (int bla = 0; bla < genome->sites.size()-1; ++bla){
      //   // variance += std::pow(site-mean, 2);
      //   variance += std::pow(genome->sites[bla]-mean, 2);
      // }
      // variance /= genome->sites.size()-1;
      // -------------
      score = 4*std::sin(variance-mean/4) + (mean*2);
      // score = 4*std::sin(4*std::sqrt(variance)-mean/4) + (mean*2);
      // score = (std::sqrt(variance) * std::sin(std::sin(mean)*variance)) + (mean*2); //max score should be 256 +/- a bit for veering off into variance land
      // (√x)*sin(sin(y)*x) + y/2
      // score = (std::tanh(variance) * std::sin(std::sin(mean)*variance)) + (mean*2);
      // score = std::sin(8*variance)*std::sqrt(variance)*std::sin(std::sin(mean)*variance) + (mean*2);
      org->dataMap.append("mean", mean); //average each mean together
      org->dataMap.append(name+"mean", mean); //record individual mean
      org->dataMap.append("variance", variance); //average each variance together
      org->dataMap.append(name+"variance", variance); //record individual variance
    }
    
    org->dataMap.append("score", score); //average each score together
		org->dataMap.append(name+"score", score); //record individual scores
    
  }
}

void
GenomicMultiBrainWorld::evaluate(std::map<std::string, std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug) {
  for (auto& org:groups[groupNamePL->get(PT)]->population){
    evaluateSolo(org, analyze, visualize, debug);
  }
}

std::unordered_map<std::string, std::unordered_set<std::string>>
GenomicMultiBrainWorld::requiredGroups() {
  // return {{groupNamePL->get(PT), {"G:" + genomeNamePL->get(PT)}}};
  std::unordered_set<std::string> blas;
	for (auto name : genomeNames) {
		blas.insert("G:" + name);
	}
	return { { groupNamePL->get(PT), blas } };
}



