#include "EncodeGenerator.h"
#include <fstream>
#include <sstream>
#include <unordered_set>
namespace kn {

static std::array<Atom, EncodeGenerator::kLengthOfEncodes> GetAtomListHelper(
    const Config &config,
    const std::pair<int, int> &jump_pair) {
  std::unordered_set<int> atom1_first_neighbors_set(
      config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborList().begin(),
      config.GetAtomList()[jump_pair.first].GetFirstNearestNeighborList().end());
  atom1_first_neighbors_set.erase(jump_pair.second);

  const std::unordered_set<int> atom1_second_neighbors_set(
      config.GetAtomList()[jump_pair.first].GetSecondNearestNeighborList().begin(),
      config.GetAtomList()[jump_pair.first].GetSecondNearestNeighborList().end());
  const std::unordered_set<int> atom1_third_neighbors_set(
      config.GetAtomList()[jump_pair.first].GetThirdNearestNeighborList().begin(),
      config.GetAtomList()[jump_pair.first].GetThirdNearestNeighborList().end());

  std::unordered_set<int> atom2_first_neighbors_set(
      config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborList().begin(),
      config.GetAtomList()[jump_pair.second].GetFirstNearestNeighborList().end());
  atom2_first_neighbors_set.erase(jump_pair.first);

  const std::unordered_set<int> atom2_second_neighbors_set(
      config.GetAtomList()[jump_pair.second].GetSecondNearestNeighborList().begin(),
      config.GetAtomList()[jump_pair.second].GetSecondNearestNeighborList().end());
  const std::unordered_set<int> atom2_third_neighbors_set(
      config.GetAtomList()[jump_pair.second].GetThirdNearestNeighborList().begin(),
      config.GetAtomList()[jump_pair.second].GetThirdNearestNeighborList().end());

  std::array<std::unordered_set<int>, EncodeGenerator::kNumOfSubEncode> sub_encode_sets;
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

  std::array<Atom, EncodeGenerator::kLengthOfEncodes> atom_list;
  int count = 0;
  for (const auto &sub_encode_set:sub_encode_sets) {
    for (const auto &index:sub_encode_set) {
      atom_list[count++] = config.GetAtomList()[index];
    }
  }
  return atom_list;
}
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
  for (const int index : first_atom.GetFirstNearestNeighborList()) {
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
static std::array<int, 59> RotateAtomsAndGetCodeHelper(
    std::array<Atom, EncodeGenerator::kLengthOfEncodes> &atom_list,
    const Matrix33 &rotation_matrix,
    const std::unordered_map<std::string, int> &type_category_hashmap) {
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
  auto first_it = atom_list.begin();
  std::array<Atom, EncodeGenerator::kLengthOfEncodes>::iterator second_it;
  for (const auto kSubCodeLength : EncodeGenerator::All_Encode_Length) {
    second_it = first_it + kSubCodeLength;
    std::sort(first_it, second_it,
              [](const Atom &first_atom, const Atom &second_atom) {
                return first_atom.GetRelativePosition() < second_atom.GetRelativePosition();
              });
    first_it = second_it;
  }

  std::array<int, EncodeGenerator::kLengthOfEncodes> encode{};
  size_t index = 0;
  for (const auto &atom : atom_list) {
    encode[index++] = type_category_hashmap.at(atom.GetType());
  }
  return encode;
}

std::vector<std::array<int, EncodeGenerator::kLengthOfEncodes>> EncodeGenerator::Encode(
    const Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, int> &type_category_hashmap) {
  auto atom_list = GetAtomListHelper(config, jump_pair);

  const auto move_distance = Vector3{0.5, 0.5, 0.5} - GetPairCenterHelper(config, jump_pair);
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
                                               GetJumpMatrixHelper(config, jump_pair),
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
std::array<int, EncodeGenerator::kLengthOfEncodes> EncodeGenerator::GetBackwardEncode(
    const std::array<int, kLengthOfEncodes> &forward_encode) {

  std::array<int, kLengthOfEncodes> backward_encode{};
  auto first_it = const_cast<int *>(forward_encode.begin());
  std::array<int, kLengthOfEncodes>::iterator second_it;
  auto insert_it = backward_encode.begin();
  for (const auto kSubCodeLength : EncodeGenerator::All_Encode_Length) {
    second_it = first_it + kSubCodeLength;
    std::reverse_copy(first_it, second_it, insert_it);
    insert_it += kSubCodeLength;
    first_it = second_it;
  }
  return backward_encode;
}

void EncodeGenerator::PrintOutEncode(
    const std::string &reference_filename,
    const std::unordered_map<std::string, int> &type_category_hashmap) {
  std::ifstream ifs(reference_filename, std::ifstream::in);
  std::ofstream ofs("encode.txt", std::ofstream::out);
  std::string buffer;
  while (ifs >> buffer) {
    if (buffer != "config") {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    int config_index, image_index, jump_pair_first, jump_pair_second;
    // config 0 end 0 pair: 248 218
    ifs >> config_index >> buffer >> image_index >> buffer >> jump_pair_first >> jump_pair_second;

    auto config = Config::ReadConfig("config" + std::to_string(config_index) + "/s/start.cfg",
                                     true);
    auto forward_encode_result =
        Encode(config, {jump_pair_first, jump_pair_second}, type_category_hashmap);
    for (const auto &image_forward_encode : forward_encode_result) {
      ofs << "config " << config_index << " end " << image_index << "  ";
      for (const auto &code : image_forward_encode) {
        ofs << code << " ";
      }
      for (const auto &code : GetBackwardEncode(image_forward_encode)) {
        ofs << code << " ";
      }
      ofs << '\n';
    }
  }
}


} // namespace kn
