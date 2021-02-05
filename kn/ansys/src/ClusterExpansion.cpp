#include "ClusterExpansion.h"
#include "Cluster.hpp"

namespace ansys::ClusterExpansion {

using Singlet_t = cfg::Cluster<1>;
using Pair_t = cfg::Cluster<2>;
using Triplet_t = cfg::Cluster<3>;

std::unordered_map<std::string, std::vector<double>> GetOneHotEncodeHashmap(
    const std::set<std::string> &type_set) {
  size_t num_singlets = type_set.size();
  std::unordered_map<std::string, std::vector<double>> encode_dict;

  size_t ct1 = 0;
  for (const auto &element : type_set) {
    std::vector<double> element_encode(num_singlets, 0);
    element_encode[ct1] = 1.0;
    encode_dict[element] = element_encode;
    ++ct1;
  }

  size_t num_pairs = num_singlets * num_singlets;
  size_t ct2 = 0;
  for (const auto &element1 : type_set) {
    for (const auto &element2 : type_set) {
      std::vector<double> element_encode(num_pairs, 0);
      element_encode[ct2] = 1.0;
      encode_dict[element1 + element2] = element_encode;
      ++ct2;
    }
  }
  return encode_dict;
}

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

// Helps to sort the atoms first symmetrically and then positionally
static bool AtomSortCompare(const cfg::Atom &lhs, const cfg::Atom &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  if (diff_x < -kEpsilon)
    return true;
  if (diff_x > kEpsilon)
    return false;

  const double diff_y_sym = std::abs(relative_position_lhs[kYDimension] - 0.5)
      - std::abs(relative_position_rhs[kYDimension] - 0.5);
  if (diff_y_sym < -kEpsilon)
    return true;
  if (diff_y_sym > kEpsilon)
    return false;

  const double diff_z_sym = std::abs(relative_position_lhs[kZDimension] - 0.5)
      - std::abs(relative_position_rhs[kZDimension] - 0.5);
  if (diff_z_sym < -kEpsilon)
    return true;
  if (diff_z_sym > kEpsilon)
    return false;
  // sort by position if they are same
  const double y_diff = relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
  if (y_diff < -kEpsilon)
    return true;
  if (y_diff > kEpsilon)
    return false;

  return relative_position_lhs[kZDimension] < relative_position_rhs[kZDimension] - kEpsilon;
}

// Helps to sort the atoms symmetrically
template<size_t DataSize>
static bool IsClusterSmallerSymmetrically(const cfg::Cluster<DataSize> &lhs,
                                          const cfg::Cluster<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    const auto &relative_position_lhs = lhs.GetAtomAt(i).GetRelativePosition();
    const auto &relative_position_rhs = rhs.GetAtomAt(i).GetRelativePosition();

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

    const double diff_z = std::abs(relative_position_lhs[kZDimension] - 0.5)
        - std::abs(relative_position_rhs[kZDimension] - 0.5);
    if (diff_z < -kEpsilon)
      return true;
    if (diff_z > kEpsilon)
      return false;
  }
  // if it reaches here, it means that the clusters are same symmetrically. Returns false.
  return false;
}

static std::vector<cfg::Atom> RotateAtomVectorAndSortHelper(
    std::vector<cfg::Atom> &&atom_list,
    const cfg::Config &reference_config,
    const std::pair<size_t, size_t> &jump_pair) {
  RotateAtomVector(atom_list, GetPairRotationMatrix(reference_config, jump_pair));
  //sort using mm2 group point
  std::sort(atom_list.begin(), atom_list.end(),
            [](const cfg::Atom &lhs, const cfg::Atom &rhs) -> bool {
              return AtomSortCompare(lhs, rhs);
            });

  size_t new_id = 0;
  for (auto &atom : atom_list) {
    atom.SetId(new_id++);
  }
  cfg::Config config(reference_config.GetBasis(), std::move(atom_list));
  //Todo only update first nearest neighbors
  config.UpdateFirstAndSecondNeighbors();
  return config.GetAtomList();
}

// Returns forward and backward sorted atom lists
static std::array<std::vector<cfg::Atom>, 2> GetSymmetricallySortedAtomVectors(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) {
  // First, second, third nearest neighbors of the jump pairs
  constexpr size_t kNumOfAtoms = 60;

  std::unordered_set<size_t>
      atom_id_hashset = GetFirstAndSecondThirdNeighborsSetOfJumpPair(config, jump_pair);

  const auto move_distance = Vector_t{0.5, 0.5, 0.5} - GetPairCenter(config, jump_pair);

  std::vector<cfg::Atom> atom_list_forward;
  atom_list_forward.reserve(kNumOfAtoms);

  Vector_t vacancy_relative_position, vacancy_cartesian_position;
  for (auto id : atom_id_hashset) {
    cfg::Atom atom = config.GetAtomList()[id];
    atom.CleanNeighborsLists();
    // move to center
    auto relative_position = atom.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);

    if (atom.GetId() == jump_pair.first) {
      vacancy_relative_position = atom.GetRelativePosition();
      vacancy_cartesian_position = atom.GetCartesianPosition();
      continue;
    }

    atom_list_forward.push_back(std::move(atom));
  }

  auto atom_list_backward(atom_list_forward);
  auto jump_atom_it_backward = std::find_if(atom_list_backward.begin(), atom_list_backward.end(),
                                            [&jump_pair](const auto &atom) {
                                              return atom.GetId() == jump_pair.second;
                                            });
  jump_atom_it_backward->SetRelativePosition(vacancy_relative_position);
  jump_atom_it_backward->SetCartesianPosition(vacancy_cartesian_position);

  return {RotateAtomVectorAndSortHelper(std::move(atom_list_forward),
                                        config,
                                        jump_pair),
          RotateAtomVectorAndSortHelper(std::move(atom_list_backward),
                                        config,
                                        {jump_pair.second, jump_pair.first})};
}

std::array<std::vector<std::string>, 2> GetForwardAndBackwardEncode(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) {
  const auto atom_vectors = GetSymmetricallySortedAtomVectors(config, jump_pair);

  // First, second, third nearest neighbors of the jump pairs
  constexpr size_t kNumOfAtoms = 60;

  std::array<std::vector<std::string>, 2> result;
  for (size_t i = 0; i < 2; ++i) {
    result[i].reserve(kNumOfAtoms);
    const auto &atom_vector = atom_vectors[i];
    std::transform(atom_vector.begin(),
                   atom_vector.end(),
                   std::back_inserter(result[i]),
                   [](const auto &atom) { return atom.GetType(); });
  }
  return result;
}

template<size_t DataSize>
static void GetAverageParametersMappingFromClusterVectorHelper(
    std::vector<cfg::Cluster<DataSize>> &&cluster_vector,
    std::vector<std::vector<std::vector<size_t>>> &cluster_mapping) {
  // sort clusters
  std::sort(cluster_vector.begin(), cluster_vector.end(),
            [](const auto &lhs, const auto &rhs) {
              return IsClusterSmallerSymmetrically(lhs, rhs);
            });
  // start to point at Cluster in the first range
  typename std::vector<cfg::Cluster<DataSize>>::const_iterator lower_it, upper_it;
  lower_it = cluster_vector.cbegin();

  do {
    upper_it = std::upper_bound(lower_it, cluster_vector.cend(),
                                *lower_it,
                                [](const auto &lhs, const auto &rhs) {
                                  return IsClusterSmallerSymmetrically(lhs, rhs);
                                });
    std::vector<std::vector<size_t>> cluster_index_vector;
    for (auto it = lower_it; it != upper_it; ++it) {
      std::vector<size_t> cluster_index;
      cluster_index.reserve(DataSize);
      for (size_t i = 0; i < DataSize; ++i) {
        cluster_index.push_back(it->GetAtomAt(i).GetId());
      }
      cluster_index_vector.push_back(cluster_index);
    }
    cluster_mapping.push_back(cluster_index_vector);

    // update to next range
    lower_it = upper_it;
  } while (upper_it != cluster_vector.cend());
}

std::vector<std::vector<std::vector<size_t>>> GetAverageClusterParametersMapping(
    const cfg::Config &config) {

  size_t vacancy_index = cfg::GetVacancyIndex(config);
  size_t neighbor_index = config.GetAtomList()[vacancy_index].GetFirstNearestNeighborsList()[0];
  const auto atom_vector = GetSymmetricallySortedAtomVectors(config, {vacancy_index,
                                                                      neighbor_index})[0];

  std::vector<std::vector<std::vector<size_t>>> cluster_mapping{};
  /// singlets
  // find all singlets
  std::vector<Singlet_t> singlet_vector;
  std::transform(atom_vector.begin(),
                 atom_vector.end(),
                 std::back_inserter(singlet_vector),
                 [](auto &&atom) {
                   return static_cast<Singlet_t>(atom);
                 });
  GetAverageParametersMappingFromClusterVectorHelper(
      std::move(singlet_vector), cluster_mapping);
  /// first nearest pairs
  // find all pairs
  std::unordered_set<Pair_t, boost::hash<Pair_t>> first_pair_set;
  for (const auto &atom1 : atom_vector) {
    for (const auto &atom2_index : atom1.GetFirstNearestNeighborsList()) {
      first_pair_set.emplace(atom1, atom_vector.at(atom2_index));
    }
  }
  GetAverageParametersMappingFromClusterVectorHelper(
      std::vector<Pair_t>(first_pair_set.begin(), first_pair_set.end()), cluster_mapping);
  /// second nearest pairs
  // find all pairs
  std::unordered_set<Pair_t, boost::hash<Pair_t>> second_pair_set;
  for (const auto &atom1 : atom_vector) {
    for (const auto &atom2_index : atom1.GetSecondNearestNeighborsList()) {
      second_pair_set.emplace(atom1, atom_vector.at(atom2_index));
    }
  }
  GetAverageParametersMappingFromClusterVectorHelper(
      std::vector<Pair_t>(second_pair_set.begin(), second_pair_set.end()), cluster_mapping);
  // /// first nearest triplets
  // // find all triplets
  // std::unordered_set<Triplet_t, boost::hash<Triplet_t>> triplets_set;
  // for (const auto &atom1 : atom_vector) {
  //   for (const auto &atom2_index : atom1.GetFirstNearestNeighborsList()) {
  //     const auto &atom2 = atom_vector[atom2_index];
  //     for (const auto &atom3_index : atom2.GetFirstNearestNeighborsList()) {
  //       if (std::find(atom1.GetFirstNearestNeighborsList().begin(),
  //                     atom1.GetFirstNearestNeighborsList().end(),
  //                     atom3_index) != atom1.GetFirstNearestNeighborsList().end()) {
  //         triplets_set.emplace(atom1, atom2, atom_vector.at(atom3_index));
  //       }
  //     }
  //   }
  // }
  // GetAverageParametersMappingFromClusterVectorHelper(
  //     std::vector<Triplet_t>(triplets_set.begin(), triplets_set.end()), cluster_mapping);
  return cluster_mapping;
}

std::vector<double> GetOneHotParametersFromMap(
    const std::vector<std::string> &encode,
    const std::unordered_map<std::string, std::vector<double>> &one_hot_encode_hashmap,
    size_t num_of_elements,
    const std::vector<std::vector<std::vector<size_t>>> &cluster_mapping) {
  constexpr size_t kEncodeListSizeThirdToDimer = 909;
  std::vector<double> res_encode;
  res_encode.reserve(kEncodeListSizeThirdToDimer);
  for (const auto &cluster_vector : cluster_mapping) {
    std::vector<double> sum_of_list(
        static_cast<size_t>(std::pow(num_of_elements, cluster_vector[0].size())), 0);
    for (const auto &cluster : cluster_vector) {
      std::string cluster_type;
      for (auto index : cluster) {
        cluster_type += encode[index];
      }
      const auto &cluster_one_hot_encode = one_hot_encode_hashmap.at(cluster_type);
      std::transform(sum_of_list.begin(), sum_of_list.end(),
                     cluster_one_hot_encode.begin(),
                     sum_of_list.begin(),
                     std::plus<>());
    }
    auto cluster_vector_size = static_cast<double>( cluster_vector.size());
    std::for_each(sum_of_list.begin(),
                  sum_of_list.end(),
                  [cluster_vector_size](auto &n) { n /= cluster_vector_size; });

    std::move(sum_of_list.begin(), sum_of_list.end(), std::back_inserter(res_encode));
  }
  return res_encode;
}
} // namespace ansys::ClusterExpansion
