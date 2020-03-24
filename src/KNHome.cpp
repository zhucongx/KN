//
// Created by Zhucong Xi on 2/11/20.
//

// #include "../include/KNHome.h"
//
// int KNHome::generateOrderedSingle(const int &i, int index,
//                                   const std::vector<int> &dupFactors,
//                                   const double &latticeConstant,
//                                   const std::string &POT,
//                                   const FCCEmbededCluster::occupInfo_256 &o256) {
//   Config cnf0;
//   cnf0.generateFCC(latticeConstant, "Al", dupFactors);
//   std::vector<std::pair<std::string, std::string>>
//       elemPairs = {{"Al", "Mg"}, {"Al", "Zn"}, {"Mg", "Zn"}};
//   const std::vector<std::pair<int, int>> &jumpPairsRef = o256.jumpPairs[i];
//   int subIndex = 0;
//   for (const auto &elemPair:elemPairs) {
//     Config cnf1 = cnf0;
//     cnf1.embedCluster(elemPair,o256,i);
//     cnf1.writeConfig(std::to_string(index)+".cfg");
//   }
//   return 0;
// }
// void KNHome::createOrdered(Config &cnfModifier,
//                            const std::vector<int> &dupFactors,
//                            const double &LC,
//                            const std::string &POT) {
//   FCCEmbededCluster::occupInfo_256 o256;
//   int index = 0;
//   for (int i = 0; i < o256.mapping.size(); ++i)
//     index = createOrderedSingle(i, index, cnfModifier, dupFactors, \
//                                 LC, POT, o256);
// }