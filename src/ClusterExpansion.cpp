#include "ClusterExpansion.h"
#include <unordered_set>
namespace kn::ClusterExpansion {
static Vector3 GetPairCenterHelper(const Config &config,
                                   const std::pair<int, int> &jump_pair) {
  Vector3 center_position;
  for (const auto kDim : All_Dimensions) {
    double first_relative = config.GetAtomList()[jump_pair.first].GetRelativePosition()[kDim];
    const double
        second_relative = config.GetAtomList()[jump_pair.second].GetRelativePosition()[kDim];

    double distance = first_relative - second_relative;
    int period = static_cast<int>(distance / 0.5);
    // make sure distance is the range (0, 0.5)
    while (period != 0) {
      first_relative -= static_cast<double>(period);
      distance = first_relative - second_relative;
      period = static_cast<int>(distance / 0.5);
    }
    center_position[kDim] = 0.5 * (first_relative + second_relative);
  }
  return center_position;
}

static Matrix33 GetJumpMatrixHelper(const Config &config,
                                    const std::pair<int, int> &jump_pair) {
  const Vector3 pair_direction = Normalize(GetRelativeDistanceVector(
      config.GetAtomList()[jump_pair.first],
      config.GetAtomList()[jump_pair.second]));
  const Atom &first_atom = config.GetAtomList()[jump_pair.first];
  Vector3 vertical_vector;
  for (const int index : first_atom.GetFirstNearestNeighborsList()) {
    const Vector3 jump_vector = GetRelativeDistanceVector(first_atom, config.GetAtomList()[index]);
    const double dot_prod = Dot(pair_direction, jump_vector);
    if (abs(dot_prod) < 1e-6) {
      double absolute_distance_square = Inner(jump_vector * config.GetBasis());
      if (absolute_distance_square < pow(Al_const::kFirstNearestNeighborsCutoff, 2)) {
        vertical_vector = Normalize(jump_vector);
        break;
      }
    }
  }
  // The third row is normalized since it is a cross product of two normalized vectors.
  // We use transposed matrix here because transpose of an orthogonal matrix equals its inverse
  return TransposeMatrix33({pair_direction, vertical_vector,
                            Cross(pair_direction, vertical_vector)});
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


static bool AreSameSymmetrically(const Atom &lhs, const Atom &rhs) {
  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();
  return abs(relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension]) < kEpsilon &&
      abs(abs(relative_position_lhs[kYDimension] - 0.5)
              - abs(relative_position_rhs[kYDimension] - 0.5)) < kEpsilon &&
      abs(abs(relative_position_lhs[kZDimension] - 0.5)
              - abs(relative_position_rhs[kZDimension] - 0.5)) < kEpsilon;
}
static bool AreSameSymmetrically(const Pair &lhs, const Pair &rhs) {
  return AreSameSymmetrically(lhs.GetAtom1(), rhs.GetAtom1()) &&
      AreSameSymmetrically(lhs.GetAtom2(), rhs.GetAtom2());
}
static bool AreSameSymmetrically(const Triplet &lhs, const Triplet &rhs) {
  return AreSameSymmetrically(lhs.GetAtom1(), rhs.GetAtom1()) &&
      AreSameSymmetrically(lhs.GetAtom2(), rhs.GetAtom2()) &&
      AreSameSymmetrically(lhs.GetAtom3(), rhs.GetAtom3());
}
// static bool AreSameSymmetrically(const Vector3 &relative_position_lhs, const Vector3 &relative_position_rhs) {
//   return abs(relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension]) < kEpsilon &&
//       abs(abs(relative_position_lhs[kYDimension] - 0.5) - abs(relative_position_rhs[kYDimension] - 0.5)) < kEpsilon &&
//       abs(abs(relative_position_lhs[kZDimension] - 0.5) - abs(relative_position_rhs[kZDimension] - 0.5)) < kEpsilon;
// }
static bool IsSmallerSymmetrically(const Atom &lhs, const Atom &rhs) {

  const auto &relative_position_lhs = lhs.GetRelativePosition();
  const auto &relative_position_rhs = rhs.GetRelativePosition();

  const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
  // const double diff_y = relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];

  if (diff_x < -kEpsilon)
    return true;
  if (diff_x > kEpsilon)
    return false;

  const double diff_y =
      abs(relative_position_lhs[kYDimension] - 0.5) - abs(relative_position_rhs[kYDimension] - 0.5);
  if (diff_y < -kEpsilon)
    return true;
  if (diff_y > kEpsilon)
    return false;

  return (abs(relative_position_lhs[kZDimension] - 0.5)
      < abs(relative_position_rhs[kZDimension] - 0.5) - kEpsilon);
}
static bool IsSmallerSymmetrically(const Pair &lhs, const Pair &rhs) {
  if (IsSmallerSymmetrically(lhs.GetAtom1(), rhs.GetAtom1()))
    return true;
  if (IsSmallerSymmetrically(rhs.GetAtom1(), lhs.GetAtom1()))
    return false;
  return IsSmallerSymmetrically(lhs.GetAtom2(), rhs.GetAtom2());
}
static bool IsSmallerSymmetrically(const Triplet &lhs, const Triplet &rhs) {
  if (IsSmallerSymmetrically(lhs.GetAtom1(), rhs.GetAtom1()))
    return true;
  if (IsSmallerSymmetrically(rhs.GetAtom1(), lhs.GetAtom1()))
    return false;
  if (IsSmallerSymmetrically(lhs.GetAtom2(), rhs.GetAtom2()))
    return true;
  if (IsSmallerSymmetrically(rhs.GetAtom2(), lhs.GetAtom2()))
    return false;
  return IsSmallerSymmetrically(lhs.GetAtom3(), rhs.GetAtom3());
}

static void RotateHelper(
    std::vector<Atom> &atom_list,
    const Matrix33 &rotation_matrix) {
  const auto move_distance_after_rotation = Vector3{0.5, 0.5, 0.5}
      - (Vector3{0.5, 0.5, 0.5} * rotation_matrix);
  for (auto &atom : atom_list) {
    auto relative_position = atom.GetRelativePosition();
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);
  }

}

static Config GetRotatedCenteredSortedConfig(const Config &config,
                                             const std::pair<int, int> &jump_pair) {
  // # First, second, third nearest neighbors of the jump pairs
  constexpr int kNumOfAtoms = 60;

  std::unordered_set<int> atom_id_set;
  atom_id_set.merge(config.GetAtomList()[jump_pair.first].GetFirstAndSecondThirdNeighborsSet());
  atom_id_set.merge(config.GetAtomList()[jump_pair.second].GetFirstAndSecondThirdNeighborsSet());

  const auto move_distance = Vector3{0.5, 0.5, 0.5} - GetPairCenterHelper(config, jump_pair);

  std::vector<Atom> atom_list;
  atom_list.reserve(kNumOfAtoms);
  for (int id : atom_id_set) {
    Atom atom = config.GetAtomList()[id];

    // move to center
    auto relative_position = atom.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);

    atom_list.push_back(std::move(atom));
  }

  RotateHelper(atom_list, GetJumpMatrixHelper(config, jump_pair));

  //sort
  std::sort(atom_list.begin(), atom_list.end(),
            [](const Atom &lhs, const Atom &rhs) {
              return IsSmallerSymmetrically(lhs, rhs);
            });

  Config config_out(config.GetBasis(), kNumOfAtoms);
  for (const Atom &atom : atom_list) {
    config_out.AppendAtomWithChangingAtomID(atom);
  }
  config_out.UpdateNeighbors();

#ifndef NDEBUG
  Config::WriteConfig(config_out, "1.cfg");
#endif

  return config_out;
}

// static void GetAverageClusterFunctionsOfSinglets(const Config &reference_config,
//                                                  );
std::vector<double> GetAverageClusterFunctions(
    const Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap) {

  Config transformed_config = GetRotatedCenteredSortedConfig(config, jump_pair);

  // Get near neighbors atom list, sorted using mm2 group point
  const auto &atom_list_reference = transformed_config.GetAtomList();

  // std::set<Atom> atom_set;
  // std::copy(atom_list.cbegin(),
  //           atom_list.cend(),
  //           std::inserter(atom_set, atom_set.end()));

  std::vector<double> average_cluster_functions_vector;

  /// singlets
  // start to point at atoms in the first range
  std::vector<Atom>::const_iterator lower_it_singlet, upper_it_singlet;
  lower_it_singlet = atom_list_reference.cbegin();
  do {
    upper_it_singlet = std::upper_bound(lower_it_singlet, atom_list_reference.cend(),
                                        *lower_it_singlet,
                                        [](const Atom &lhs, const Atom &rhs) {
                                          return IsSmallerSymmetrically(lhs, rhs);
                                        });
    // add singlets
    average_cluster_functions_vector.push_back(
        std::accumulate(lower_it_singlet, upper_it_singlet, 0.0,
                        [&type_category_hashmap = std::as_const(type_category_hashmap)](
                            double current_sum, const Atom &rhs) {
                          return current_sum
                              + type_category_hashmap.at(rhs.GetType());
                        }) / std::distance(lower_it_singlet, upper_it_singlet)
    );
    // update to next range
    lower_it_singlet = upper_it_singlet;
  } while (upper_it_singlet != atom_list_reference.cend());

  /// pairs
  // find all pairs
  std::array<std::unordered_set<Pair, boost::hash<Pair>>, 2> pairs_sets;
  // std::unordered_set<Pair, boost::hash<Pair> > pairs_set;

  for (const auto &atom1 : atom_list_reference) {
    for (const auto &atom2_index : atom1.GetFirstNearestNeighborsList()) {
      pairs_sets[0].emplace(atom1, atom_list_reference[atom2_index]);
      // pairs_set.emplace(atom1, atom_list_reference[atom2_index]);
    }
    for (const auto &atom2_index : atom1.GetSecondNearestNeighborsList()) {
      pairs_sets[1].emplace(atom1, atom_list_reference[atom2_index]);
    }
    // for (const auto &atom2_index : atom1.GetThirdNearestNeighborsList()) {
    //   pairs_sets[2].emplace(atom1, atom_list_reference[atom2_index]);
    // }
  }
  std::array<std::vector<Pair>, 2>
      neighbors_vectors{std::vector<Pair>(pairs_sets[0].begin(), pairs_sets[0].end()),
                        std::vector<Pair>(pairs_sets[1].begin(), pairs_sets[1].end()),
      // std::vector<Pair>(pairs_sets[2].begin(), pairs_sets[2].end())
  };
  // std::vector<Pair> pairs_vector(pairs_set.begin(), pairs_set.end());
  for (auto &pairs_vector : neighbors_vectors) {
    // sort pairs
    std::sort(pairs_vector.begin(), pairs_vector.end(),
              [](const Pair &lhs, const Pair &rhs) {
                return IsSmallerSymmetrically(lhs, rhs);
              });

    // start to point at atoms in the first range
    std::vector<Pair>::const_iterator lower_it_pair, upper_it_pair;
    lower_it_pair = pairs_vector.cbegin();
    do {
      upper_it_pair = std::upper_bound(lower_it_pair, pairs_vector.cend(),
                                       *lower_it_pair,
                                       [](const Pair &lhs, const Pair &rhs) {
                                         return IsSmallerSymmetrically(lhs, rhs);
                                       });
      // add pairs
      average_cluster_functions_vector.push_back(
          std::accumulate(lower_it_pair, upper_it_pair, 0.0,
                          [&type_category_hashmap = std::as_const(type_category_hashmap)](
                              double current_sum, const Pair &rhs) {
                            return current_sum
                                + type_category_hashmap.at(rhs.GetAtom1().GetType())
                                    * type_category_hashmap.at(rhs.GetAtom2().GetType());
                          }) / std::distance(lower_it_pair, upper_it_pair));
      // update to next range
      lower_it_pair = upper_it_pair;
    } while (upper_it_pair != pairs_vector.cend());
  }

#ifndef NDEBUG
  // for (const auto &a : neighbors_vectors[0]) {
  //   std::cout << a.GetAtom1().GetRelativePosition() - Vector3{0.5, 0.5, 0.5}
  //             << "  ||  " << a.GetAtom2().GetRelativePosition() - Vector3{0.5, 0.5, 0.5}
  //             << '\n';
  // }
#endif
  /// triplets quadruplets quintuplets
  // find all triplets
  std::unordered_set<Triplet, boost::hash<Triplet>> triplets_set;
  for (const auto &atom1 : atom_list_reference) {
    for (const auto &atom2_index : atom1.GetFirstNearestNeighborsList()) {
      const auto &atom2 = atom_list_reference[atom2_index];
      for (const auto &atom3_index : atom2.GetFirstNearestNeighborsList()) {
        if (std::find(atom1.GetFirstNearestNeighborsList().begin(),
                      atom1.GetFirstNearestNeighborsList().end(),
                      atom3_index) != atom1.GetFirstNearestNeighborsList().end()) {
          triplets_set.emplace(atom1, atom2, atom_list_reference[atom3_index]);
        }
      }
    }
  }
  std::vector<Triplet> triplets_vector(triplets_set.begin(), triplets_set.end());
  // sort triplets
  std::sort(triplets_vector.begin(), triplets_vector.end(),
            [](const Triplet &lhs, const Triplet &rhs) {
              return IsSmallerSymmetrically(lhs, rhs);
            });

  // start to point at atoms in the first range
  std::vector<Triplet>::const_iterator lower_it_triplet, upper_it_triplet;
  lower_it_triplet = triplets_vector.cbegin();
  do {
    upper_it_triplet = std::upper_bound(lower_it_triplet, triplets_vector.cend(),
                                        *lower_it_triplet,
                                        [](const Triplet &lhs, const Triplet &rhs) {
                                          return IsSmallerSymmetrically(lhs, rhs);
                                        });
    // add triplets
    average_cluster_functions_vector.push_back(
        std::accumulate(lower_it_triplet, upper_it_triplet, 0.0,
                        [&type_category_hashmap = std::as_const(type_category_hashmap)](
                            double current_sum, const Triplet &rhs) {
                          return current_sum
                              + type_category_hashmap.at(rhs.GetAtom1().GetType())
                                  * type_category_hashmap.at(rhs.GetAtom2().GetType())
                                  * type_category_hashmap.at(rhs.GetAtom3().GetType());
                        }) / std::distance(lower_it_triplet, upper_it_triplet));
    // update to next range
    lower_it_triplet = upper_it_triplet;
  } while (upper_it_triplet != triplets_vector.cend());

  return average_cluster_functions_vector;
}

std::vector<double> GetAverageClusterFunctionsBack(
    const Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap) {
  auto back_config = config;
  back_config.AtomsJump(jump_pair.first, jump_pair.second);
  return GetAverageClusterFunctions(back_config, jump_pair, type_category_hashmap);

}
} // namespace kn::ClusterExpansion