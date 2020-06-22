#include "ConfigUtility.h"
#include <unordered_set>
#include <ConfigIO.h>
namespace kn {

std::map<Bond, int> ConfigUtility::CountAllBonds(Config &config) {
  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  auto atom_list = config.GetAtomList();
  for (const auto &atom : atom_list) {
    type1 = atom.type_;
    for (const auto &atom2_id : atom.first_nearest_neighbor_list_) {
      bonds_count_map[Bond{type1, atom_list[atom2_id].type_}]++;
    }
  }
  for (auto &bond_count : bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
Vector3 GetPairCenterHelper(const Config &config,
                            const std::pair<int, int> &jump_pair) {
  Vector3 center_position;
  for (const auto kDim : All_Dimensions) {
    double first_relative = config.GetAtomList()[jump_pair.first].relative_position_[kDim];
    const double second_relative = config.GetAtomList()[jump_pair.second].relative_position_[kDim];

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
Matrix33 GetJumpMatrixHelper(const Config &config,
                             const std::pair<int, int> &jump_pair) {
  const Vector3 pair_direction = Normalize(AtomUtility::GetRelativeDistanceVector(
      config.GetAtomList()[jump_pair.first],
      config.GetAtomList()[jump_pair.second]));
  const Atom &first_atom = config.GetAtomList()[jump_pair.first];
  Vector3 vertical_vector;
  for (const int index : first_atom.first_nearest_neighbor_list_) {
    const Vector3 jump_vector = AtomUtility::GetRelativeDistanceVector(first_atom,
                                                                       config.GetAtomList()[index]);
    const double dot_prod = Dot(pair_direction, jump_vector);
    if (abs(dot_prod) < 1e-6) {
      double absolute_distance_square = Inner(jump_vector * config.GetBasis());
      if (absolute_distance_square < pow(Al_const::kFirstNearestNeighborCutoff, 2)) {
        vertical_vector = Normalize(jump_vector);
        break;
      }
    }
  }

  return {-pair_direction, -vertical_vector, Normalize(Cross(pair_direction, vertical_vector))};
}
std::vector<std::string> RotateAtomsAndGetCodeHelper(const std::string &jump_type,
                                                     std::vector<Atom> &atom_list,
                                                     const Matrix33 &rotation_matrix) {
  const auto move_distance_after_rotation = Vector3{0.5, 0.5, 0.5}
      - (Vector3{0.5, 0.5, 0.5} * rotation_matrix);
  for (auto &atom : atom_list) {
    // rotate
    atom.relative_position_ = atom.relative_position_ * rotation_matrix;
    // move to new center
    atom.relative_position_ += move_distance_after_rotation;
    atom.relative_position_ -= ElementFloor(atom.relative_position_);
  }
  std::sort(atom_list.begin(), atom_list.end(),
            [](const Atom &first_atom, const Atom &second_atom) {
              return first_atom.relative_position_ < second_atom.relative_position_;
            });
  std::vector<std::string> string_codes;
  string_codes.reserve(Al_const::kLengthOfEncodes);
  string_codes.push_back(jump_type);
  for (const auto &atom : atom_list) {
    string_codes.push_back(atom.type_);
  }

  return string_codes;
}
std::vector<std::vector<std::string>> ConfigUtility::Encode(const Config &config,
                                                            const std::pair<int, int> &jump_pair) {
  std::vector<std::vector<std::string>> result;

  std::unordered_set<int> id_set;
  for (const auto &atom:{config.GetAtomList()[jump_pair.first],
                         config.GetAtomList()[jump_pair.second]}) {
    for (const int neighbor_index : atom.first_nearest_neighbor_list_)
      id_set.insert(neighbor_index);
    for (const int neighbor_index : atom.second_nearest_neighbor_list_)
      id_set.insert(neighbor_index);
  }
  std::vector<int> temporary_id;
  std::copy_if(id_set.begin(), id_set.end(), std::back_inserter(temporary_id),
               [jump_pair](const int &index) {
                 return index != jump_pair.first && index != jump_pair.second;
               });

  std::vector<Atom> atom_list;
  for (int index : temporary_id) {
    atom_list.push_back(config.GetAtomList()[index]);
  }

  std::sort(atom_list.begin(), atom_list.end(),
            [](const Atom &first_atom, const Atom &second_atom) {
              return first_atom.relative_position_ < second_atom.relative_position_;
            });

  const auto move_distance = Vector3{0.5, 0.5, 0.5} - GetPairCenterHelper(config, jump_pair);
  for (auto &atom : atom_list) {
    // move to center
    atom.relative_position_ += move_distance;
    atom.relative_position_ -= ElementFloor(atom.relative_position_);
  }

  const auto &jump_type = config.GetAtomList()[jump_pair.second].type_;

  // First Rotation
  result.push_back(RotateAtomsAndGetCodeHelper(jump_type, atom_list,
                                               InverseMatrix33(GetJumpMatrixHelper(config,
                                                                                   jump_pair))));

  // 2-fold rotation
  result.push_back(RotateAtomsAndGetCodeHelper(jump_type, atom_list,
                                               {{{1.0, 0.0, 0.0},
                                                 {0.0, -1.0, 0.0},
                                                 {0.0, 0.0, -1.0}}}));

  // mirror y
  result.push_back(RotateAtomsAndGetCodeHelper(jump_type, atom_list,
                                               {{{1.0, 0.0, 0.0},
                                                 {0.0, -1.0, 0.0},
                                                 {0.0, 0.0, 1.0}}}));
  // mirror z
  result.push_back(RotateAtomsAndGetCodeHelper(jump_type, atom_list,
                                               {{{1.0, 0.0, 0.0},
                                                 {0.0, 1.0, 0.0},
                                                 {0.0, 0.0, -1.0}}}));

  return result;
}

} // namespace kn
