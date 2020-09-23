#ifndef KN_KN_KMC_INCLUDE_MPICLUSTERS_H_
#define KN_KN_KMC_INCLUDE_MPICLUSTERS_H_

#include "MpiIterator.h"
namespace kn {

class MpiClusters : public MpiIterator {
 public:
  MpiClusters(long long int initial_number,
              long long int increment_number,
              long long int finial_number,
              std::string solvent_element,
              int smallest_cluster_criteria,
              int solvent_bond_criteria);
  void IterateToRun() override;
 private:
  std::string solvent_element_;
  int smallest_cluster_criteria_{};
  int solvent_bond_criteria_{};
};

}// namespace kn


#endif //KN_KN_KMC_INCLUDE_MPICLUSTERS_H_
