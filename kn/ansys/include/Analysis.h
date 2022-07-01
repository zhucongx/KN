#ifndef KN_KN_KMC_INCLUDE_MPICLUSTERS_H_
#define KN_KN_KMC_INCLUDE_MPICLUSTERS_H_

#include "MpiIterator.h"
#include "ClustersFinder.h"

namespace ansys {

class Analysis{
  public:
    Analysis(unsigned long long int initial_number,
             unsigned long long int increment_number,
             std::string solvent_element,
             size_t smallest_cluster_criteria,
             size_t solvent_bond_criteria);
    ~Analysis();
    void SerialRunCluster() const;
    void SerialRunWarrenCowley() const;

  private:
    const unsigned long long initial_number_;
    const unsigned long long increment_number_;
    unsigned long long final_number_;
    std::string solvent_element_;
    size_t smallest_cluster_criteria_;
    size_t solvent_bond_criteria_;
    std::unordered_map<unsigned long long, double> filename_time_hashset_{};
};

}// namespace ansys


#endif //KN_KN_KMC_INCLUDE_MPICLUSTERS_H_
