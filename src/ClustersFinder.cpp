#include "ClustersFinder.h"
#include <queue>
#include <unordered_map>
#include <utility>
#include "ConfigIO.h"

namespace kn {
ClustersFinder::ClustersFinder(std::string cfg_filename,
                               std::string solvent_atom_type,
                               int smallest_cluster_criteria,
                               int solvent_bond_criteria,
                               double first_nearest_neighbors_distance,
                               double second_nearest_neighbors_distance)
  : cfg_filename_(std::move(cfg_filename)),
    solvent_element_(std::move(solvent_atom_type)),
    smallest_cluster_criteria_(smallest_cluster_criteria),
    solvent_bond_criteria_(solvent_bond_criteria) {
  ReadFileAndUpdateNeighbor(first_nearest_neighbors_distance,
                            second_nearest_neighbors_distance);
}

ClustersFinder::ClusterElementNumMap ClustersFinder::FindClustersAndOutput() {
  auto cluster_to_atom_map = FindAtomListOfClusters();

  Config config_out;
  config_out.SetScale(config_.GetScale());
  config_out.SetBasis(config_.GetBasis());
  std::vector<std::map<std::string, int>> num_atom_in_clusters_set;
  auto atoms_list_reference = config_.GetAtomList();
  for (auto &atom_list : cluster_to_atom_map) {
    // initialize map with all the element, because some cluster may not have all types of element
    std::map<std::string, int> num_atom_in_one_cluster;
    for (const auto &element : element_set_) {
      num_atom_in_one_cluster[element] = 0;
    }

    for (const auto &atom_index : atom_list) {
      num_atom_in_one_cluster[atoms_list_reference[atom_index].type_]++;
      config_out.AppendAtom(atoms_list_reference[atom_index]);
    }

    num_atom_in_clusters_set.push_back(std::move(num_atom_in_one_cluster));
  }
  auto const pos = cfg_filename_.find_last_of('.');
  const std::string output_name_suffix = cfg_filename_.substr(pos + 1);
  std::string output_name = cfg_filename_.substr(0, pos);
  output_name += "_cluster.";
  output_name += output_name_suffix;
  ConfigIO::WriteConfig(config_out, output_name, false);
  return num_atom_in_clusters_set;
}

void ClustersFinder::PrintLog(const std::string &filename,
                              const ClusterElementNumMap &found_data) {
  std::ofstream ofs("clusters_info.txt", std::ofstream::out | std::ofstream::app);
  ofs << '#' << filename << '\n';
  for (const auto &cluster : found_data) {
    for (const auto &[key, count] : cluster) {
      ofs << count << ' ';
    }
    ofs << '\n';
  }
}

void ClustersFinder::ReadFileAndUpdateNeighbor(double first_nearest_neighbors_distance,
                                               double second_nearest_neighbors_distance) {
  config_ = ConfigIO::ReadConfig(cfg_filename_);
  config_.UpdateNeighbors(first_nearest_neighbors_distance, second_nearest_neighbors_distance);
  for (const auto &[element_type, index_vector] : config_.GetElementListMap()) {
    if (element_type == "X")
      continue;
    element_set_.insert(element_type);
  }
}

std::unordered_set<int> ClustersFinder::FindSoluteAtomsHelper() const {
  std::unordered_set<int> solute_atoms_hashset;
  for (const auto &[element_type, index_vector] : config_.GetElementListMap()) {
    if (element_type == solvent_element_ || element_type == "X")
      continue;
    std::copy(index_vector.begin(),
              index_vector.end(),
              std::inserter(solute_atoms_hashset, solute_atoms_hashset.end()));
  }

  return solute_atoms_hashset;
}

std::vector<std::vector<int>> ClustersFinder::FindAtomListOfClustersBFSHelper(
    std::unordered_set<int> unvisited_atoms_id_set) const {
  auto atoms_list_reference = config_.GetAtomList();
  std::vector<std::vector<int>> cluster_atom_list;
  std::queue<int> visit_id_queue;
  int atom_id;

  std::unordered_set<int>::iterator it;
  while (!unvisited_atoms_id_set.empty()) {
    // Find next element
    it = unvisited_atoms_id_set.begin();
    visit_id_queue.push(*it);
    unvisited_atoms_id_set.erase(it);

    std::vector<int> atom_list_of_one_cluster;
    while (!visit_id_queue.empty()) {
      atom_id = visit_id_queue.front();
      visit_id_queue.pop();

      atom_list_of_one_cluster.push_back(atom_id);
      for (const auto &neighbor_id : atoms_list_reference[atom_id].first_nearest_neighbor_list_) {
        it = unvisited_atoms_id_set.find(neighbor_id);
        if (it != unvisited_atoms_id_set.end()) {
          visit_id_queue.push(*it);
          unvisited_atoms_id_set.erase(it);
        }
      }
    }
    cluster_atom_list.push_back(atom_list_of_one_cluster);
  }
  return cluster_atom_list;
}

std::vector<std::vector<int>> ClustersFinder::FindAtomListOfClusters() const {
  auto cluster_atom_list = FindAtomListOfClustersBFSHelper(FindSoluteAtomsHelper());

  // remove small clusters
  auto it = cluster_atom_list.begin();
  while (it != cluster_atom_list.end()) {
    if (static_cast<int>(it->size()) <= smallest_cluster_criteria_) {
      it = cluster_atom_list.erase(it);
    } else {
      ++it;
    }
  }

  // add solvent neighbors
  auto atoms_list_reference = config_.GetAtomList();
  for (auto &atom_list : cluster_atom_list) {
    std::unordered_map<int, int> neighbor_bond_count;
    for (const auto &atom_index : atom_list) {
      for (auto neighbor_id : atoms_list_reference[atom_index].first_nearest_neighbor_list_) {
        if (atoms_list_reference[neighbor_id].type_ == solvent_element_)
          neighbor_bond_count[neighbor_id]++;
      }
    }
    for (auto [neighbor_id, bond_count] : neighbor_bond_count) {
      if (bond_count >= solvent_bond_criteria_)
        atom_list.push_back(neighbor_id);
    }
  }

  // remove duplicate outer layer
  std::unordered_set<int> unvisited_atoms_id_set;
  for (const auto &kIndexVector : cluster_atom_list) {
    std::copy(kIndexVector.begin(),
              kIndexVector.end(),
              std::inserter(unvisited_atoms_id_set, unvisited_atoms_id_set.end()));
  }
  auto cluster_atom_list_after_removing = FindAtomListOfClustersBFSHelper(unvisited_atoms_id_set);

  return cluster_atom_list_after_removing;
}

} // namespace kn
