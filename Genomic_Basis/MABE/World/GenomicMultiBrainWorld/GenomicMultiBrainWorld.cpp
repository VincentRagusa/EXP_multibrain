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

std::shared_ptr<ParameterLink<bool>>
GenomicMultiBrainWorld::randomizeGenomesPL = Parameters::register_parameter( "WORLD_GMB-randomizeGenomes", (bool) 0, "seed genomes with all 0s or random values?");

std::shared_ptr<ParameterLink<double>>
GenomicMultiBrainWorld::APeriodPL = Parameters::register_parameter( "WORLD_GMB-APeriod", (double) 1.0, "how long A runs before switching to B");

std::shared_ptr<ParameterLink<double>>
GenomicMultiBrainWorld::BPeriodPL = Parameters::register_parameter( "WORLD_GMB-BPeriod", (double) 1.0, "how long B runs before switching back to A");

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::AFunctionPL = Parameters::register_parameter("WORLD_GMB-AFunction", (std::string) "VECT[0,0]", "(MTree)");

std::shared_ptr<ParameterLink<std::string>>
GenomicMultiBrainWorld::BFunctionPL = Parameters::register_parameter("WORLD_GMB-BFunction", (std::string) "VECT[0,0]", "(MTree)");

std::shared_ptr<ParameterLink<bool>>
GenomicMultiBrainWorld::seedGenomeBetweenAandBPL = Parameters::register_parameter( "WORLD_GMB-seedGenomeBetweenAandB", (bool) 0, "seed genomes with all 'APeriod' if true. The value of 'APeriod' will be down-cast as an int before it is inserted into the genome.");

std::shared_ptr<ParameterLink<bool>>
GenomicMultiBrainWorld::recordFirstValleyCrossPL = Parameters::register_parameter( "WORLD_GMB-recordFirstValleyCross", (bool) 0, "records the first generation when any task");


GenomicMultiBrainWorld::GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_):
  AbstractWorld(PT_) {
  convertCSVListToVector(genomeNamePL->get(PT), genomeNames);
  AFunctionMT = stringToMTree(AFunctionPL->get(PT));
  BFunctionMT = stringToMTree(BFunctionPL->get(PT));
  La = APeriodPL->get(PT);
  Lb = BPeriodPL->get(PT);
  randomizeGenome = randomizeGenomesPL->get(PT);
  seedGenomeBetweenAandB = seedGenomeBetweenAandBPL->get(PT);
  recordFirstValleyCross = recordFirstValleyCrossPL->get(PT);
  
  if (seedGenomeBetweenAandBPL){
    std::cout << "WORLD: WARNING!\n\tThe value of 'APeriodPL' will be down-cast as an int! This can cause missalignment between mutation size and position in the fitness landscape." << std::endl;
  }

  if (randomizeGenome && seedGenomeBetweenAandB){
    std::cout << "WORLD: ERROR!\n\tYou cannot set both 'randomizeGenome' and 'seedGenomeBetweenAandB' parameters to 1 at the same time!" << std::endl;
    exit(1);
  }

  // columns to be added to ave file
  popFileColumns.clear();
  popFileColumns.push_back("score");
  for (auto name : genomeNames) { 
		popFileColumns.push_back(name + "score");
	}

  popFileColumns.push_back("mean");
  // popFileColumns.push_back("variance");
  for (auto name : genomeNames) { 
    popFileColumns.push_back(name + "mean");
    // popFileColumns.push_back(name + "variance");
  }
  
  
}

template<typename T>
auto
arithmetic_mean(std::vector<T> data){
  return std::accumulate(data.begin(), data.end(), 0.0) / (data.size());
}

template<typename T>
auto
arithmetic_variance(std::vector<T> data, double MA){
  auto VA = 0.0;
  for (auto val:data){
    VA += std::pow(val - MA, 2);
  }
  VA /= data.size();
  return VA;
}

std::vector<double>
GenomicMultiBrainWorld::C(double x){
  double c_a = 0;
  double c_b = 0;
  while (x >= La + Lb){
    ++c_a; ++c_b;
    x -= La + Lb;
  }
  if (x >= La){
    ++c_a;
    x -= La;
  }
  return {c_a, c_b, x};
}

double
GenomicMultiBrainWorld::F(double x){
  auto c = C(x);
  auto c_a = c[0];
  auto c_b = c[1];
  auto d = c[2];


  std::vector<std::vector<double>> vec_a, vec_b;

  if (c_a == c_b){
    vec_a = {{c_a*La + d}};
    vec_b = {{c_b*Lb}};
  }
  else{
    vec_a = {{c_a*La}};
    vec_b = {{c_b*Lb + d}};
  }

  auto f_a = AFunctionMT->eval(vec_a)[0];
  auto f_b = BFunctionMT->eval(vec_b)[0];

  return f_a + f_b;
}

void
GenomicMultiBrainWorld::evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug) {
  for (auto name : genomeNames) {

    if (randomizeGenome && Global::update == 0){
      org->genomes[name]->fillRandom();  
    }
    
    auto genome = std::dynamic_pointer_cast<CircularGenome<int>>(org->genomes[name]);
    if (genome == nullptr){
      std::cout << "WORLD: ERROR!\n\tThis world requires that you run MABE with circular genome and that the 'GENOME-sitesType' parameter is set to 'int'." << std::endl;
      exit(1);
    }

    if (seedGenomeBetweenAandB && Global::update == 0){
      for (auto& site : genome->sites){
        site = (int) La;
      }
    }

    auto mean = arithmetic_mean(genome->sites);
    // auto variance = arithmetic_variance(genome->sites, mean);
    auto score = F(mean );
    org->dataMap.append("mean", mean); //average each mean together
    org->dataMap.append(name+"mean", mean); //record individual mean

    if (name == "A::"){
      org->dataMap.append("score", mean);
      org->dataMap.append(name+"score", mean);
    }
    else {
      org->dataMap.append("score", score);
      org->dataMap.append(name+"score", score);
    }
    
    
    // org->dataMap.append("variance", variance); //average each variance together
    // org->dataMap.append(name+"variance", variance); //record individual variance

    if (recordFirstValleyCross){
      if (mean >= 2*La + Lb){
        FileManager::writeToFile("valley_cross_time.csv", "-1," + std::to_string(Global::update) + "," + name, "bottom_time, top_time, trait_name");
        local_finished = true;
      }

      

      if (mean > La + Lb){ // if crossed 
        std::vector<std::string>::iterator it;
        it = std::find(passed.begin(), passed.end(), name);
        name_had_1_passed[name] |= true; // make sure we know someone with this trait is crossed so we don't reset the passed flag.

        if (it == passed.end()){ //if not already recorded, record.
          FileManager::writeToFile("valley_cross_time.csv", std::to_string(Global::update) + ",-1," + name, "bottom_time, top_time, trait_name");
          passed.push_back(name); //disable further saving
        }
      }
    }
  }
}

void
GenomicMultiBrainWorld::evaluate(std::map<std::string, std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug) {
  // at the start of every generation (call to evaluate) reset the one-passed detector.
  if (recordFirstValleyCross){
    for (auto name:genomeNames){
      name_had_1_passed[name] = false;
    }
  }

  for (auto& org:groups[groupNamePL->get(PT)]->population){
    evaluateSolo(org, analyze, visualize, debug);
  }

  if (local_finished){ // if in the previous loop someone set local_finish, tell the archavist.
    groups[groupNamePL->get(PT)]->archivist->finished_ = true;
  }
  
  if (recordFirstValleyCross){
    for (auto name:genomeNames){
      if (! name_had_1_passed[name]){ //if this trait had no pop members passed the valley
        std::vector<std::string>::iterator it;
        it = std::find(passed.begin(), passed.end(), name);

        if (it != passed.end()){ // and they are recorded as having passed previously
          passed.erase(it); //lol... erase it!
        }
      }
    }
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



