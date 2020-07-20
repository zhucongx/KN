#include "ClustersFinder.h"
#include <queue>
#include <unordered_map>
#include <utility>

namespace kn {
ClustersFinder::ClustersFinder(std::string cfg_filename,
                               std::string solvent_atom_type,
                               int smallest_cluster_criteria,
                               int solvent_bond_criteria)
    : cfg_filename_(std::move(cfg_filename)),
      solvent_element_(std::move(solvent_atom_type)),
      smallest_cluster_criteria_(smallest_cluster_criteria),
      solvent_bond_criteria_(solvent_bond_criteria) {
  ReadFileAndUpdateNeighbor();
}

ClustersFinder::ClusterElementNumMap ClustersFinder::FindClustersAndOutput() {
  auto cluster_to_atom_map = FindAtomListOfClusters();

  Config config_out(config_.GetBasis(), config_.GetNumAtoms());
  std::vector<std::map<std::string, int>> num_atom_in_clusters_set;
  for (auto &atom_list : cluster_to_atom_map) {
    // initialize map with all the element, because some cluster may not have all types of element
    std::map<std::string, int> num_atom_in_one_cluster;
    for (const auto &element : element_set_) {
      num_atom_in_one_cluster[element] = 0;
    }

    for (const auto &atom_index : atom_list) {
      num_atom_in_one_cluster[config_.GetAtomList()[atom_index].GetType()]++;
      config_out.AppendAtom(config_.GetAtomList()[atom_index]);
    }

    num_atom_in_clusters_set.push_back(std::move(num_atom_in_one_cluster));
  }
  auto output_name(cfg_filename_);
  auto const pos = output_name.find_last_of('.');
  output_name.insert(pos,"_cluster" );
  Config::WriteConfig(config_out, output_name, false);
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

void ClustersFinder::ReadFileAndUpdateNeighbor() {
  config_ = Config::ReadConfig(cfg_filename_, true);
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
      for (const auto &neighbor_id : config_.GetAtomList()[atom_id].GetFirstNearestNeighborList()) {
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
  for (auto &atom_list : cluster_atom_list) {
    std::unordered_map<int, int> neighbor_bond_count;
    for (const auto &atom_index : atom_list) {
      for (auto neighbor_id : config_.GetAtomList()[atom_index].GetFirstNearestNeighborList()) {
        if (config_.GetAtomList()[neighbor_id].GetType() == solvent_element_)
          neighbor_bond_count[neighbor_id]++;
      }
    }
    for (auto[neighbor_id, bond_count] : neighbor_bond_count) {
      if (bond_count >= solvent_bond_criteria_)
        atom_list.push_back(neighbor_id);
    }
  }

  // remove duplicate outer layer
  std::unordered_set<int> unvisited_atoms_id_set;
  for (const auto &singe_cluster_vector : cluster_atom_list) {
    std::copy(singe_cluster_vector.begin(),
              singe_cluster_vector.end(),
              std::inserter(unvisited_atoms_id_set, unvisited_atoms_id_set.end()));
  }
  auto cluster_atom_list_after_removing = FindAtomListOfClustersBFSHelper(unvisited_atoms_id_set);

  return cluster_atom_list_after_removing;
}

} // namespace kn
