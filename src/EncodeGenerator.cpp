#include "EncodeGenerator.h"
#include <fstream>
#include <sstream>
#include <unordered_set>
#include "ConfigIO.h"
namespace kn {
EncodeGenerator::EncodeGenerator(std::string reference_filename) :
    reference_filename_(std::move(reference_filename)) {
}

Vector3 GetPairCenterHelper(const Config &config,
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
Matrix33 GetJumpMatrixHelper(const Config &config,
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
      if (absolute_distance_square < pow(Al_const::kFirstNearestNeighborCutoff, 2)) {
        vertical_vector = Normalize(jump_vector);
        break;
      }
    }
  }

  return {-pair_direction, -vertical_vector, Normalize(Cross(pair_direction, vertical_vector))};
}
std::array<std::string, Al_const::kLengthOfEncodes> RotateAtomsAndGetCodeHelper(
    const std::string &jump_type,
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
  std::sort(atom_list.begin(), atom_list.end(),
            [](const Atom &first_atom, const Atom &second_atom) {
              return first_atom.GetRelativePosition() < second_atom.GetRelativePosition();
            });
  std::array<std::string, Al_const::kLengthOfEncodes> string_codes;
  string_codes[0] = jump_type;
  size_t index = 0;
  for (const auto &atom : atom_list) {
    string_codes[++index] = atom.GetType();
  }

  return string_codes;
}
std::vector<std::array<std::string, Al_const::kLengthOfEncodes>> EncodeGenerator::Encode(
    const Config &config,
    const std::pair<int, int> &jump_pair) {

  // put jump type and environment atoms' index into a set to avoid duplicate
  std::unordered_set<int> id_set;
  for (const auto &atom:{config.GetAtomList()[jump_pair.first],
                         config.GetAtomList()[jump_pair.second]}) {
    for (const int neighbor_index : atom.GetFirstNearestNeighborList())
      id_set.insert(neighbor_index);
    for (const int neighbor_index : atom.GetNearNeighborList())
      id_set.insert(neighbor_index);
  }
  // convert set to vector and don't keep the jump type
  std::vector<int> temporary_id;
  temporary_id.reserve(Al_const::kLengthOfEncodes);
  std::copy_if(id_set.begin(), id_set.end(), std::back_inserter(temporary_id),
               [jump_pair](const int &index) {
                 return index != jump_pair.first && index != jump_pair.second;
               });

  std::vector<Atom> atom_list;
  atom_list.reserve(temporary_id.size());
  for (int index : temporary_id) {
    atom_list.push_back(config.GetAtomList()[index]);
  }

  std::sort(atom_list.begin(), atom_list.end(),
            [](const Atom &first_atom, const Atom &second_atom) {
              return first_atom.GetRelativePosition() < second_atom.GetRelativePosition();
            });

  const auto move_distance = Vector3{0.5, 0.5, 0.5} - GetPairCenterHelper(config, jump_pair);
  for (auto &atom : atom_list) {
    // move to center
    auto relative_position = atom.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);

    atom.SetRelativePosition(relative_position);
  }

  const auto &jump_type = config.GetAtomList()[jump_pair.second].GetType();
  std::vector<std::array<std::string, Al_const::kLengthOfEncodes>> result;
  result.reserve(4);
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
void EncodeGenerator::PrintOutEncode() {
  std::ifstream ifs(reference_filename_, std::ifstream::in);
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

    auto config = ConfigIO::ReadConfig("config" + std::to_string(config_index) + "/s/start.cfg",
                                       true);
    auto encode_result = Encode(config, {jump_pair_first, jump_pair_second});
    for (const auto &image_encode : encode_result) {
      ofs << "config " << config_index << " end " << image_index << "   ";
      for (const auto &code : image_encode) {
        ofs << code << "   ";
      }
      ofs << '\n';
    }
  }
}

} // namespace kn
