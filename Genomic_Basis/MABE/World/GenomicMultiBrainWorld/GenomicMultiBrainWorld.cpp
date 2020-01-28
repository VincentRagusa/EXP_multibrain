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


GenomicMultiBrainWorld::GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_):
  AbstractWorld(PT_) {
  convertCSVListToVector(genomeNamePL->get(PT), genomeNames);
  AFunctionMT = stringToMTree(AFunctionPL->get(PT));
  BFunctionMT = stringToMTree(BFunctionPL->get(PT));
  La = APeriodPL->get(PT);
  Lb = BPeriodPL->get(PT);
  randomizeGenome = randomizeGenomesPL->get(PT);


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
      std::cout << "GARBO set your genome site type to int ya fool (not really i'm the fool)" << std::endl;
      exit(1);
    }

    auto mean = arithmetic_mean(genome->sites);
    // auto variance = arithmetic_variance(genome->sites, mean);
    auto score = F(mean );
    
    org->dataMap.append("score", score);
    org->dataMap.append(name+"score", score);
    org->dataMap.append("mean", mean); //average each mean together
    org->dataMap.append(name+"mean", mean); //record individual mean
    // org->dataMap.append("variance", variance); //average each variance together
    // org->dataMap.append(name+"variance", variance); //record individual variance
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



