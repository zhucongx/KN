#include "ClusterExpansion.h"
#include "Cluster.hpp"

namespace ansys::ClusterExpansion {

using Singlet_t = cfg::Cluster<1>;
using Pair_t = cfg::Cluster<2>;
using Triplet_t = cfg::Cluster<3>;
// mm2 group point:
// mirror perpendicular to xz plane
// {1 0 0
//  0 -1 0
//  0 0 1}
// mirror perpendicular to xy plane
// {1 0 0
//  0 1 0
//  0 0 -1}
// 2 fold rotation along x axis
// {1 0 0
//  0 -1 0
//  0 0 -1}
/// When we treat this vector, we don't care the sign of y and z, and y and z should be
/// indistinguishable. If sum and diff of abs of y and abs of z are same, these atoms should
/// be considered same positions
// Helps to sort the atoms symmetrically
bool IsAtomSmallerSymmetrically(const cfg::Atom &lhs, const cfg::Atom &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];

  if (diff_x < -kEpsilon)
    return true;
  if (diff_x > kEpsilon)
    return false;

  const double diff_y = std::abs(relative_position_lhs[kYDimension] - 0.5)
      - std::abs(relative_position_rhs[kYDimension] - 0.5);
  if (diff_y < -kEpsilon)
    return true;
  if (diff_y > kEpsilon)
    return false;

  return (std::abs(relative_position_lhs[kZDimension] - 0.5)
      < std::abs(relative_position_rhs[kZDimension] - 0.5) - kEpsilon);
}

// Helps to sort the atoms symmetrically
template<size_t DataSize>
static bool IsClusterSmallerSymmetrically(const cfg::Cluster<DataSize> &lhs,
                                          const cfg::Cluster<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {

    if (IsAtomSmallerSymmetrically(lhs.GetAtomAt(i), rhs.GetAtomAt(i)))
      return true;
    if (IsAtomSmallerSymmetrically(rhs.GetAtomAt(i), lhs.GetAtomAt(i)))
      return false;
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}

// Get the neighboring updated sorted atoms vector
static std::vector<cfg::Atom> GetSymmetricallySortedAtomVector(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) {
  // First, second, third nearest neighbors of the jump pairs
  constexpr size_t kNumOfAtoms = 60;

  std::unordered_set<size_t> atom_id_set;
  atom_id_set.merge(config.GetAtomList()[jump_pair.first].GetFirstAndSecondThirdNeighborsSet());
  atom_id_set.merge(config.GetAtomList()[jump_pair.second].GetFirstAndSecondThirdNeighborsSet());
  atom_id_set.erase(jump_pair.first);

  const auto move_distance = Vector_t{0.5, 0.5, 0.5} - GetPairCenter(config, jump_pair);

  std::vector<cfg::Atom> atom_list;
  atom_list.reserve(kNumOfAtoms);
  for (auto id : atom_id_set) {
    cfg::Atom atom = config.GetAtomList()[id];
    atom.CleanNeighborsLists();
    // move to center
    auto relative_position = atom.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);

    atom_list.push_back(std::move(atom));
  }
  RotateAtomVector(atom_list, GetPairRotationMatrix(config, jump_pair));

  //sort using mm2 group point
  std::sort(atom_list.begin(), atom_list.end(),
            [](const auto &lhs, const auto &rhs) {
              return IsAtomSmallerSymmetrically(lhs, rhs);
            });

  size_t new_id = 0;
  for (auto &atom : atom_list) {
    atom.SetId(new_id++);
  }

  cfg::Config config_out(config.GetBasis(), std::move(atom_list));
  config_out.UpdateNeighbors();

  return config_out.GetAtomList();
}

// Update the average cluster functions vector from the cluster vector
template<size_t DataSize>
static void GetAverageParametersFromClusterVectorHelper(
    std::vector<cfg::Cluster<DataSize>> &&cluster_vector,
    const std::unordered_map<std::string, double> &type_category_hashmap,
    std::vector<double> &average_cluster_functions_vector) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) {
              return IsClusterSmallerSymmetrically(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<cfg::Cluster<DataSize>>::const_iterator lower_it_pair, upper_it_pair;
  lower_it_pair = cluster_vector.cbegin();
  do {
    upper_it_pair = std::upper_bound(lower_it_pair, cluster_vector.cend(),
                                     *lower_it_pair,
                                     [](const auto &lhs, const auto &rhs) {
                                       return IsClusterSmallerSymmetrically(lhs, rhs);
                                     });
    // add parameters
    average_cluster_functions_vector.push_back(
        std::accumulate(lower_it_pair, upper_it_pair, 0.0,
                        [&type_category_hashmap = std::as_const(type_category_hashmap)](
                            double current_sum, const cfg::Cluster<DataSize> &rhs) {
                          double cumulative_product = 1.0;
                          for (size_t i = 0; i < DataSize; ++i) {
                            cumulative_product *=
                                type_category_hashmap.at(rhs.GetAtomAt(i).GetType());
                          }
                          return current_sum + cumulative_product;
                        }) / static_cast<double>(std::distance(lower_it_pair, upper_it_pair))
    );

    // update to next range
    lower_it_pair = upper_it_pair;
  } while (upper_it_pair != cluster_vector.cend());
}

std::vector<double> GetAverageClusterParameters(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap) {

  const auto atom_vector = GetSymmetricallySortedAtomVector(config, jump_pair);
  std::vector<double> average_cluster_functions_vector;

  /// singlets
  // find all singlets
  std::vector<Singlet_t> singlet_vector;
  std::transform(atom_vector.begin(),
                 atom_vector.end(),
                 std::back_inserter(singlet_vector),
                 [](auto &&atom) {
                   return static_cast<Singlet_t>(atom);
                 });
  GetAverageParametersFromClusterVectorHelper(std::move(singlet_vector),
                                              type_category_hashmap,
                                              average_cluster_functions_vector);

  /// first nearest pairs
  // find all pairs
  std::unordered_set<Pair_t, boost::hash<Pair_t>> first_pair_set;
  for (const auto &atom1 : atom_vector) {
    for (const auto &atom2_index : atom1.GetFirstNearestNeighborsList()) {
      first_pair_set.emplace(atom1, atom_vector.at(atom2_index));
    }
  }
  GetAverageParametersFromClusterVectorHelper(
      std::vector<Pair_t>(first_pair_set.begin(), first_pair_set.end()),
      type_category_hashmap,
      average_cluster_functions_vector);

  /// second nearest pairs
  // find all pairs
  std::unordered_set<Pair_t, boost::hash<Pair_t>> second_pair_set;
  for (const auto &atom1 : atom_vector) {
    for (const auto &atom2_index : atom1.GetSecondNearestNeighborsList()) {
      second_pair_set.emplace(atom1, atom_vector.at(atom2_index));
    }
  }
  GetAverageParametersFromClusterVectorHelper(
      std::vector<Pair_t>(second_pair_set.begin(), second_pair_set.end()),
      type_category_hashmap,
      average_cluster_functions_vector);

  /// first nearest triplets
  // find all triplets
  std::unordered_set<Triplet_t, boost::hash<Triplet_t>> triplets_set;
  for (const auto &atom1 : atom_vector) {
    for (const auto &atom2_index : atom1.GetFirstNearestNeighborsList()) {
      const auto &atom2 = atom_vector[atom2_index];
      for (const auto &atom3_index : atom2.GetFirstNearestNeighborsList()) {
        if (std::find(atom1.GetFirstNearestNeighborsList().begin(),
                      atom1.GetFirstNearestNeighborsList().end(),
                      atom3_index) != atom1.GetFirstNearestNeighborsList().end()) {
          triplets_set.emplace(atom1, atom2, atom_vector.at(atom3_index));
        }
      }
    }
  }
  GetAverageParametersFromClusterVectorHelper(
      std::vector<Triplet_t>(triplets_set.begin(), triplets_set.end()),
      type_category_hashmap,
      average_cluster_functions_vector);

  return average_cluster_functions_vector;
}

std::vector<double> GetAverageClusterParametersBack(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap) {
  auto back_config = config;
  cfg::AtomsJump(back_config, jump_pair.first, jump_pair.second);
  return GetAverageClusterParameters(back_config, jump_pair, type_category_hashmap);
}
} // namespace ansys::ClusterExpansion
