#ifndef KN_KN_KMC_INCLUDE_MPICLUSTERS_H_
#define KN_KN_KMC_INCLUDE_MPICLUSTERS_H_

#include "MpiIterator.h"
namespace kn {

class MpiClusters : public MpiIterator {
  public:
    MpiClusters(unsigned long long int initial_number,
                unsigned long long int increment_number,
                unsigned long long int finial_number,
                std::string solvent_element,
                size_t smallest_cluster_criteria,
                size_t solvent_bond_criteria);
    virtual ~MpiClusters();
    void IterateToRun() override;
  private:
    std::string solvent_element_;
    size_t smallest_cluster_criteria_{};
    size_t solvent_bond_criteria_{};
};

}// namespace kn


#endif //KN_KN_KMC_INCLUDE_MPICLUSTERS_H_