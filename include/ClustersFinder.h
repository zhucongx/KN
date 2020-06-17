#ifndef KN_INCLUDE_CLUSTERSFINDER_H_
#define KN_INCLUDE_CLUSTERSFINDER_H_
#include <vector>
#include <set>
#include <unordered_set>

#include "FCCConfig.h"

namespace kn {
class ClustersFinder {
  public:
    using ClusterElementNumMap = std::vector<std::map<std::string, int>>;
    ClustersFinder(std::string cfg_filename,
                   std::string solvent_atom_type,
                   int smallest_cluster_criteria,
                   int solvent_bond_criteria);
    // Return a 2D array where values of each row representing the number of atoms of different
    // element in one cluster
    ClusterElementNumMap FindClustersAndOutput();

    static void PrintLog(const std::string &filename, const ClusterElementNumMap &found_data);
  private:
    void ReadFileAndUpdateNeighbor();
    [[nodiscard]] std::unordered_set<int> FindSoluteAtomsHelper() const;
    [[nodiscard]] std::vector<std::vector<int>> FindAtomListOfClustersBFSHelper(
        std::unordered_set<int> unvisited_atoms_id_set) const;
    [[nodiscard]] std::vector<std::vector<int>> FindAtomListOfClusters() const;

    std::string cfg_filename_;
    Config config_;
    std::string solvent_element_;
    std::set<std::string> element_set_;
    int smallest_cluster_criteria_{};
    int solvent_bond_criteria_{};
};
} // namespace kn


#endif //KN_INCLUDE_CLUSTERSFINDER_H_
