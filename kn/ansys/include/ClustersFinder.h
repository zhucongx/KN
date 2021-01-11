#ifndef KN_KN_KMC_INCLUDE_CLUSTERSFINDER_H_
#define KN_KN_KMC_INCLUDE_CLUSTERSFINDER_H_
#include <vector>
#include <set>
#include <unordered_set>

#include "Config.h"

namespace ansys {
class ClustersFinder {
  public:
    using ClusterElementNumMap = std::vector<std::map<std::string, size_t>>;
    ClustersFinder(std::string cfg_filename,
                   std::string solvent_atom_type,
                   size_t smallest_cluster_criteria,
                   size_t solvent_bond_criteria);

    ClusterElementNumMap FindClustersAndOutput();

    // static void PrintLog(const unsigned long long int &file_index,
    //                      double time,
    //                      const ClusterElementNumMap &found_data);
  private:
    void UpdateElementSet();
    [[nodiscard]] std::unordered_set<size_t> FindSoluteAtomsHelper() const;
    [[nodiscard]] std::vector<std::vector<size_t>> FindAtomListOfClustersBFSHelper(
        std::unordered_set<size_t> unvisited_atoms_id_set) const;

    // Return a 2D array where values of each row representing the number of atoms of different
    // element in one cluster
    [[nodiscard]] std::vector<std::vector<size_t>> FindAtomListOfClusters() const;

    std::string cfg_filename_;
    cfg::Config config_;
    std::string solvent_element_;
    std::set<std::string> element_set_{};
    size_t smallest_cluster_criteria_{};
    size_t solvent_bond_criteria_{};
};
} // namespace ansys


#endif //KN_KN_KMC_INCLUDE_CLUSTERSFINDER_H_
