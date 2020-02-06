//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once

#include "../AbstractWorld.h"

class GenomicMultiBrainWorld : public AbstractWorld {

public:

  static std::shared_ptr<ParameterLink<std::string>> groupNamePL;
  static std::shared_ptr<ParameterLink<std::string>> genomeNamePL;
  static std::shared_ptr<ParameterLink<bool>> randomizeGenomesPL;
  static std::shared_ptr<ParameterLink<double>> APeriodPL;
  static std::shared_ptr<ParameterLink<double>> BPeriodPL;
  static std::shared_ptr<ParameterLink<std::string>> AFunctionPL;
  static std::shared_ptr<ParameterLink<std::string>> BFunctionPL;
  static std::shared_ptr<ParameterLink<bool>> seedGenomeBetweenAandBPL;
  static std::shared_ptr<ParameterLink<bool>> recordFirstValleyCrossPL;

  std::shared_ptr<Abstract_MTree> AFunctionMT, BFunctionMT;
  double La, Lb;
  bool randomizeGenome, seedGenomeBetweenAandB, recordFirstValleyCross;
  bool local_finished = false;

  double
  F(double x );

  std::vector<double>
  C(double x );
 
  std::vector<int> site_symbol_counts = std::vector<int>(256, 0);
  std::vector<std::string> genomeNames;

  GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_ = nullptr);

  virtual ~GenomicMultiBrainWorld() = default;

  void
  evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug);

  void
  evaluate(std::map<std::string,std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug);

  virtual std::unordered_map<std::string, std::unordered_set<std::string>>
  requiredGroups() override;
};

