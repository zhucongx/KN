#include "ClusterFinder.h"
#include <queue>
#include <unordered_map>
#include <utility>
namespace box {
ClusterFinder::ClusterFinder(std::string cfg_file_name,
                             std::string solvent_atom_type,
                             int smallest_cluster_criteria,
                             int solvent_bond_criteria,
                             double first_nearest_neighbors_distance,
                             double second_nearest_neighbors_distance)
    : cfg_file_name_(std::move(cfg_file_name)),
      solvent_atom_type_(std::move(solvent_atom_type)),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria) {
  ReadFileAndUpdateNeighbor(first_nearest_neighbors_distance,
                            second_nearest_neighbors_distance);
}

void ClusterFinder::ReadFileAndUpdateNeighbor(double first_nearest_neighbors_distance,
                                              double second_nearest_neighbors_distance) {
  config_.ReadConfig(cfg_file_name_);
  config_.UpdateNeighbors(first_nearest_neighbors_distance, second_nearest_neighbors_distance);
}

std::unordered_set<int> ClusterFinder::FindSoluteAtomsHelper() {
  std::unordered_set<int> solute_atoms_hashset;
  for (const auto &[element_type, index_vector] : config_.GetElementListMap()) {
    if (element_type == solvent_atom_type_ || element_type == "X")
      continue;
    std::copy(index_vector.begin(),
              index_vector.end(),
              std::inserter(solute_atoms_hashset, solute_atoms_hashset.end()));
  }

  return solute_atoms_hashset;
}
std::vector<std::vector<int>> ClusterFinder::FindAtomListOfClustersBFSHelper(
    std::unordered_set<int> unvisited_atoms_id_set) {

  std::vector<std::vector<int>> cluster_atom_list;
  std::queue<int> visit_id_queue;
  int atom_id;

  std::unordered_set<int>::iterator it;
  while (!unvisited_atoms_id_set.empty()) {
    // Find next element
    it = unvisited_atoms_id_set.begin();
    visit_id_queue.push(*it);
    unvisited_atoms_id_set.erase(it);

    std::vector<int> atom_list;
    while (!visit_id_queue.empty()) {
      atom_id = visit_id_queue.front();
      visit_id_queue.pop();

      atom_list.push_back(atom_id);
      for (const auto &neighbor_id:config_.GetAtom(atom_id).first_nearest_neighbor_list_) {
        it = unvisited_atoms_id_set.find(neighbor_id);
        if (it != unvisited_atoms_id_set.end()) {
          visit_id_queue.push(*it);
          unvisited_atoms_id_set.erase(it);
        }
      }
    }
    cluster_atom_list.push_back(atom_list);
  }
  return cluster_atom_list;
}
std::unordered_set<int> ClusterFinder::ConvertClusterAtomListToHashSetHelper(
    const std::vector<std::vector<int>> &cluster_atom_list) {
  std::unordered_set<int> unvisited_atoms_id_set;
  for (const auto &index_vector:cluster_atom_list) {
    std::copy(index_vector.begin(),
              index_vector.end(),
              std::inserter(unvisited_atoms_id_set, unvisited_atoms_id_set.end()));
  }

  return unvisited_atoms_id_set;
}
std::vector<std::vector<int>> ClusterFinder::FindAtomListOfClusters() {
  auto cluster_atom_list_all = FindAtomListOfClustersBFSHelper(FindSoluteAtomsHelper());

  // remove small clusters
  std::vector<std::vector<int>> cluster_atom_list;
  for (auto &&atom_list : cluster_atom_list_all) {
    if (atom_list.size() > smallest_cluster_criteria_)
      cluster_atom_list.push_back(std::move(atom_list));
  }

  // add solvent neighbors
  for (auto &atom_list : cluster_atom_list) {
    std::unordered_map<int, int> neighbor_bond_count;
    for (const auto &atom_index:atom_list) {
      for (auto neighbor_id:config_.GetAtom(atom_index).first_nearest_neighbor_list_) {
        if (config_.GetAtom(neighbor_id).GetType() == solvent_atom_type_)
          neighbor_bond_count[neighbor_id]++;
      }
    }
    for (auto[neighbor_id, bond_count]:neighbor_bond_count) {
      if (bond_count >= solvent_bond_criteria_)
        atom_list.push_back(neighbor_id);
    }
  }

  // remove duplicate outer layer
  auto cluster_atom_list_after_removing =
      FindAtomListOfClustersBFSHelper(ConvertClusterAtomListToHashSetHelper(cluster_atom_list));

  return cluster_atom_list_after_removing;
}

std::vector<int> ClusterFinder::FindClustersAndOutput() {
  auto cluster_to_atom_map = FindAtomListOfClusters();
  Config config_out;
  config_out.SetScale(config_.GetScale());
  config_out.SetBasis(config_.GetBasis());
  for (auto &atom_list:cluster_to_atom_map) {
    for (const auto &atom_index:atom_list) {
      config_out.AppendAtom(config_.GetAtom(atom_index));
    }
  }
  auto const pos = cfg_file_name_.find_last_of('.');
  const std::string output_name_suffix = cfg_file_name_.substr(pos + 1);
  std::string output_name = cfg_file_name_.substr(0, pos);
  output_name += "_cluster.";
  output_name += output_name_suffix;

  config_out.WriteConfig(output_name);
  std::vector<int> num_atom_in_clusters_list(config_out.GetElementListMap().size());
  for (const auto &[atom, atom_list]:config_out.GetElementListMap()) {
    num_atom_in_clusters_list.push_back(atom_list.size());
  }

  return num_atom_in_clusters_list;
}

}// namespace box
