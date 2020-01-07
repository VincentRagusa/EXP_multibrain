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
  static std::shared_ptr<ParameterLink<int>> modePL;
  static std::shared_ptr<ParameterLink<int>> numberOfOutputsPL;
  static std::shared_ptr<ParameterLink<int>> evaluationsPerGenerationPL;

  static std::shared_ptr<ParameterLink<std::string>> groupNamePL;
  static std::shared_ptr<ParameterLink<std::string>> brainNamePL;

  GenomicMultiBrainWorld(std::shared_ptr<ParametersTable> PT_ = nullptr);

  virtual ~GenomicMultiBrainWorld() = default;

  void
  evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug);

  void
  evaluate(std::map<std::string, std::shared_ptr<Group>> &groups, int analyze, int visualize, int debug);

  virtual std::unordered_map<std::string, std::unordered_set<std::string>>
  requiredGroups() override;
};

