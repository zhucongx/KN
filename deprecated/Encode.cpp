#include "Encode.h"

namespace kn::Encode {
static std::array<cfg::Atom, kLengthOfEncodes> GetAtomListHelper(
    const cfg::Config &config,
    const std::pair<int, int> &jump_pair) {
  std::unordered_set<int> atom1_first_neighbors_set(
      config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborsList().begin(),
      config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborsList().end());
  atom1_first_neighbors_set.erase(jump_pair.second);

  const std::unordered_set<int> atom1_second_neighbors_set(
      config.GetAtomList()[jump_pair.first].GetSecondNearestNeighborsList().begin(),
      config.GetAtomList()[jump_pair.first].GetSecondNearestNeighborsList().end());
  const std::unordered_set<int> atom1_third_neighbors_set(
      config.GetAtomList()[jump_pair.first].GetThirdNearestNeighborsList().begin(),
      config.GetAtomList()[jump_pair.first].GetThirdNearestNeighborsList().end());

  std::unordered_set<int> atom2_first_neighbors_set(
      config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborsList().begin(),
      config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborsList().end());
  atom2_first_neighbors_set.erase(jump_pair.first);

  const std::unordered_set<int> atom2_second_neighbors_set(
      config.GetAtomList()[jump_pair.second].GetSecondNearestNeighborsList().begin(),
      config.GetAtomList()[jump_pair.second].GetSecondNearestNeighborsList().end());
  const std::unordered_set<int> atom2_third_neighbors_set(
      config.GetAtomList()[jump_pair.second].GetThirdNearestNeighborsList().begin(),
      config.GetAtomList()[jump_pair.second].GetThirdNearestNeighborsList().end());

  std::array<std::unordered_set<int>, kNumOfSubEncode> sub_encode_sets;
  sub_encode_sets[0] = {jump_pair.second};
  for (const int index:atom1_first_neighbors_set) {
    if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
      sub_encode_sets[1].insert(index);
      continue;
    }
    if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
      sub_encode_sets[2].insert(index);
      continue;
    }
    if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
      sub_encode_sets[3].insert(index);
      continue;
    }
    sub_encode_sets[4].insert(index);
  }
  for (const int index:atom2_first_neighbors_set) {
    if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
      continue;
    }
    if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
      sub_encode_sets[2].insert(index);
      continue;
    }
    if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
      sub_encode_sets[3].insert(index);
      continue;
    }
    sub_encode_sets[4].insert(index);
  }

  for (const int index:atom1_second_neighbors_set) {
    if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
      continue;
    }
    if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
      sub_encode_sets[5].insert(index);
      continue;
    }
    if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
      sub_encode_sets[6].insert(index);
      continue;
    }
    sub_encode_sets[7].insert(index);
  }
  for (const int index:atom2_second_neighbors_set) {
    if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
      continue;
    }
    if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
      continue;
    }
    if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
      sub_encode_sets[6].insert(index);
      continue;
    }
    sub_encode_sets[7].insert(index);
  }

  for (const int index:atom1_third_neighbors_set) {
    if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
      continue;
    }
    if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
      continue;
    }
    if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
      sub_encode_sets[8].insert(index);
      continue;
    }
    sub_encode_sets[9].insert(index);
  }
  for (const int index:atom2_third_neighbors_set) {
    if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
      continue;
    }
    if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
      continue;
    }
    if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
      continue;
    }
    sub_encode_sets[9].insert(index);
  }

  std::array<cfg::Atom, kLengthOfEncodes> atom_list;
  int count = 0;
  for (const auto &sub_encode_set:sub_encode_sets) {
    for (const auto &index:sub_encode_set) {
      atom_list[count++] = config.GetAtomList()[index];
    }
  }
  return atom_list;
}

template<size_t DataSize>
static void RotateAndSort(
    std::array<cfg::Atom, DataSize> &atom_list,
    const Matrix_t &rotation_matrix) {
  const auto move_distance_after_rotation = Vector_t{0.5, 0.5, 0.5}
      - (Vector_t{0.5, 0.5, 0.5} * rotation_matrix);
  for (auto &atom : atom_list) {
    auto relative_position = atom.GetRelativePosition();
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);
  }
  //sort
  auto first_it = atom_list.begin();
  typename std::array<cfg::Atom, DataSize>::iterator second_it;
  for (auto kSubCodeLength:All_Encode_Length) {
    second_it = first_it + kSubCodeLength;
    std::sort(first_it, second_it,
              [](const cfg::Atom &first_atom, const cfg::Atom &second_atom) {
                return first_atom.GetRelativePosition() < second_atom.GetRelativePosition();
              });
    first_it = second_it;
    // stop in case that using a part of neighbors
    if (first_it >= atom_list.end())
      break;
  }
}

static std::array<int, kLengthOfEncodes> RotateAtomsAndGetCodeHelper(
    std::array<cfg::Atom, kLengthOfEncodes> &atom_list,
    const Matrix_t &rotation_matrix,
    const std::unordered_map<std::string, int> &type_category_hashmap) {

  RotateAndSort(atom_list, rotation_matrix);
  // encode
  std::array<int, kLengthOfEncodes> encode{};
  size_t index = 0;
  for (const auto &atom : atom_list) {
    encode[index++] = type_category_hashmap.at(atom.GetType());
  }
  return encode;
}

std::vector<std::array<int, kLengthOfEncodes>> GetEncode(
    const cfg::Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, int> &type_category_hashmap) {
  auto atom_list = GetAtomListHelper(config, jump_pair);

  const auto move_distance = Vector_t{0.5, 0.5, 0.5} - GetPairCenter(config, jump_pair);
  for (auto &atom : atom_list) {
    // move to center
    auto relative_position = atom.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);
  }

  std::vector<std::array<int, kLengthOfEncodes>> result;
  result.reserve(4);
  // First Rotation
  result.push_back(RotateAtomsAndGetCodeHelper(atom_list,
                                               GetPairRotationMatrix(config, jump_pair),
                                               type_category_hashmap));

  // 2-fold rotation
  result.push_back(RotateAtomsAndGetCodeHelper(atom_list,
                                               {{{1.0, 0.0, 0.0},
                                                 {0.0, -1.0, 0.0},
                                                 {0.0, 0.0, -1.0}}}, type_category_hashmap));

  // mirror y
  result.push_back(RotateAtomsAndGetCodeHelper(atom_list,
                                               {{{1.0, 0.0, 0.0},
                                                 {0.0, -1.0, 0.0},
                                                 {0.0, 0.0, 1.0}}}, type_category_hashmap));

  // mirror z
  result.push_back(RotateAtomsAndGetCodeHelper(atom_list,
                                               {{{1.0, 0.0, 0.0},
                                                 {0.0, 1.0, 0.0},
                                                 {0.0, 0.0, -1.0}}}, type_category_hashmap));

  return result;
}
std::array<int, kLengthOfEncodes> GetBackwardEncode(
    const std::array<int, kLengthOfEncodes> &forward_encode) {

  std::array<int, kLengthOfEncodes> backward_encode{};
  auto first_it = const_cast<int *>(forward_encode.begin());
  std::array<int, kLengthOfEncodes>::iterator second_it;
  auto insert_it = backward_encode.begin();
  for (const auto kSubCodeLength : All_Encode_Length) {
    second_it = first_it + kSubCodeLength;
    std::reverse_copy(first_it, second_it, insert_it);
    insert_it += kSubCodeLength;
    first_it = second_it;
  }
  return backward_encode;
}

// static std::array<Atom, kLengthOfFirstNeighbors> GetFirstAtomListHelper(
//     const Config &config,
//     const std::pair<int, int> &jump_pair) {
//   std::unordered_set<int> atom1_first_neighbors_set(
//       config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborsList().end());
//   atom1_first_neighbors_set.erase(jump_pair.second);
//
//   const std::unordered_set<int> atom1_second_neighbors_set(
//       config.GetAtomList()[jump_pair.first].GetSecondNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.first].GetSecondNearestNeighborsList().end());
//   const std::unordered_set<int> atom1_third_neighbors_set(
//       config.GetAtomList()[jump_pair.first].GetThirdNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.first].GetThirdNearestNeighborsList().end());
//
//   std::unordered_set<int> atom2_first_neighbors_set(
//       config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborsList().end());
//   atom2_first_neighbors_set.erase(jump_pair.first);
//
//   const std::unordered_set<int> atom2_second_neighbors_set(
//       config.GetAtomList()[jump_pair.second].GetSecondNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.second].GetSecondNearestNeighborsList().end());
//   const std::unordered_set<int> atom2_third_neighbors_set(
//       config.GetAtomList()[jump_pair.second].GetThirdNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.second].GetThirdNearestNeighborsList().end());
//
//   std::array<std::unordered_set<int>, 5> sub_encode_sets;
//   sub_encode_sets[0] = {jump_pair.second};
//   for (const int index:atom1_first_neighbors_set) {
//     if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
//       sub_encode_sets[1].insert(index);
//       continue;
//     }
//     if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
//       sub_encode_sets[2].insert(index);
//       continue;
//     }
//     if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
//       sub_encode_sets[3].insert(index);
//       continue;
//     }
//     sub_encode_sets[4].insert(index);
//   }
//   for (const int index:atom2_first_neighbors_set) {
//     if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
//       continue;
//     }
//     if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
//       sub_encode_sets[2].insert(index);
//       continue;
//     }
//     if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
//       sub_encode_sets[3].insert(index);
//       continue;
//     }
//     sub_encode_sets[4].insert(index);
//   }
//
//   std::array<Atom, kLengthOfFirstNeighbors> atom_list;
//   int count = 0;
//   for (const auto &sub_encode_set:sub_encode_sets) {
//     for (const auto &index:sub_encode_set) {
//       atom_list[count++] = config.GetAtomList()[index];
//     }
//   }
//   return atom_list;
// }

//
// std::pair<std::unordered_map<Bond, int>, std::unordered_map<Bond, int>> GetBondAroundPair(
//     const Config &config,
//     const std::pair<int, int> &jump_pair) {
//
//   std::unordered_set<int> atom1_first_neighbors_set(
//       config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborsList().end());
//   atom1_first_neighbors_set.erase(jump_pair.second);
//
//   std::unordered_set<int> atom2_first_neighbors_set(
//       config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborsList().begin(),
//       config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborsList().end());
//   atom2_first_neighbors_set.erase(jump_pair.first);
//
//   std::unordered_map<Bond, int> bonds_around_first, bonds_around_second;
//
//   const auto &atom_list = config.GetAtomList();
//   for (const int index:atom1_first_neighbors_set) {
//     if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
//       continue;
//     }
//     bonds_around_first[{atom_list[jump_pair.second].GetType(), atom_list[index].GetType()}]++;
//   }
//   for (const int index:atom2_first_neighbors_set) {
//     if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
//       continue;
//     }
//     bonds_around_second[{atom_list[jump_pair.second].GetType(), atom_list[index].GetType()}]++;;
//   }
//   return std::make_pair(bonds_around_first, bonds_around_second);
// }

// std::unordered_map<std::string, double> GetFirstNearestEnvironment(
//     const kn::Config &config,
//     const std::pair<int, int> &jump_pair) {
//   auto atom_list = GetFirstAtomListHelper(config, jump_pair);
//
//   const auto move_distance = Vector_t{0.5, 0.5, 0.5} - GetPairCenterHelper(config, jump_pair);
//   for (auto &atom : atom_list) {
//     // move to center
//     auto relative_position = atom.GetRelativePosition();
//     relative_position += move_distance;
//     relative_position -= ElementFloor(relative_position);
//
//     atom.SetRelativePosition(relative_position);
//   }
//   const auto rotation_matrix = GetJumpMatrixHelper(config, jump_pair);
//   RotateAndSort(atom_list, rotation_matrix);
//
//   // std::unordered_map<std::string, Vector_t> element_mean_position_hashmap;
//   // decompose the force along x direction and radius direction
//   // for (auto length:First_Neighbors_Encode_length) {
//   //   for (int i = 0; i < length / 2; i++) {
//   //
//   //   }
//   // }
//
//   return std::unordered_map<std::string, double>();
// }
} // namespace kn::Encode
