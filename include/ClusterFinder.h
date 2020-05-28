#ifndef KN_INCLUDE_CLUSTERFINDER_H_
#define KN_INCLUDE_CLUSTERFINDER_H_
#include <unordered_set>
#include <vector>

#include "FCCConfig.h"
namespace box {

class ClusterFinder {
 public:
  ClusterFinder(std::string cfg_file_name,
                std::string solvent_atom_type,
                int smallest_cluster_criteria,
                int solvent_bond_criteria,
                double first_nearest_neighbors_distance = Al_const::kFirstNearestNeighborCutoff,
                double second_nearest_neighbors_distance = Al_const::kSecondNearestNeighborsCutoff);

  void FindClustersAndOutput();
 private:
  void ReadFileAndUpdateNeighbor(double first_nearest_neighbors_distance,
                                 double second_nearest_neighbors_distance);
  std::unordered_set<int> FindSoluteAtomsHelper();
  std::vector<std::vector<int>> FindAtomListOfClustersBFSHelper(std::unordered_set<int> unvisited_atoms_id_set);
  std::unordered_set<int> ConvertClusterAtomListToHashSetHelper(const std::vector<std::vector<int>> &cluster_atom_list);

  std::vector<std::vector<int>> FindAtomListOfClusters();


  std::string cfg_file_name_;
  Config config_;
  std::string solvent_atom_type_;
  int smallest_cluster_criteria_{};
  int solvent_bond_criteria_{};

};

}// namespace box


#endif //KN_INCLUDE_CLUSTERFINDER_H_
