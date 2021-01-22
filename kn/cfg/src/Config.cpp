#include"Config.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <boost/range/combine.hpp>
namespace cfg {

Config::Config() = default;
Config::Config(const Matrix_t &basis, size_t atom_size) : basis_(basis) {
  if (atom_size)
    atom_list_.reserve(atom_size);
}
Config::Config(const Matrix_t &basis, std::vector<Atom> atom_list)
    : basis_(basis), atom_list_(std::move(atom_list)) {}
size_t Config::GetNumAtoms() const {
  return atom_list_.size();
}
const Matrix_t &Config::GetBasis() const {
  return basis_;
}
const std::vector<Atom> &Config::GetAtomList() const {
  return atom_list_;
}
std::map<std::string, std::vector<size_t>> Config::GetElementListMap() const {
  std::map<std::string, std::vector<size_t>> element_list_map;
  for (const auto &atom : atom_list_) {
    element_list_map[atom.GetType()].push_back(atom.GetId());
  }
  return element_list_map;
}
std::set<std::string> Config::GetTypeSet() const {
  std::set<std::string> type_set;
  for (const auto &atom : atom_list_) {
    type_set.insert(atom.GetType());
  }
  return type_set;
}
void Config::ConvertRelativeToCartesian() {
  for (auto &atom : atom_list_) {
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}
void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix(basis_);
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetCartesianPosition() * inverse_basis);
  }
}
void Config::ScaleWith(double scale) {
  basis_ *= scale;
  ConvertRelativeToCartesian();
}
void Config::WrapAtomRelative() {
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(
        atom.GetRelativePosition() - ElementFloor(atom.GetRelativePosition()));
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}
void Config::WrapAtomCartesian() {
  auto inverse_basis = InverseMatrix(basis_);
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetCartesianPosition() * inverse_basis);
    atom.SetRelativePosition(
        atom.GetRelativePosition() - ElementFloor(atom.GetRelativePosition()));
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}
void Config::MoveRelativeDistance(const Vector_t &distance_vector) {
  for (auto &atom : atom_list_) {
    atom.SetRelativePosition(atom.GetRelativePosition() + distance_vector);
    atom.SetRelativePosition(
        atom.GetRelativePosition() - ElementFloor(atom.GetRelativePosition()));
    atom.SetCartesianPosition(atom.GetRelativePosition() * basis_);
  }
}
void Config::MoveOneAtomRelativeDistance(size_t index,
                                         const Vector_t &distance_vector) {
  atom_list_[index].SetRelativePosition(
      atom_list_[index].GetRelativePosition() + distance_vector);
  atom_list_[index].SetRelativePosition(atom_list_[index].GetRelativePosition()
                                            - ElementFloor(atom_list_[index].GetRelativePosition()));
  atom_list_[index].SetCartesianPosition(atom_list_[index].GetRelativePosition() * basis_);
}
void Config::Perturb(std::mt19937_64 &generator) {
  constexpr double kMean = 0;
  constexpr double kStandardDeviation = 0.15;
  constexpr double kPerturbCutOff = 0.3;
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);

  for (auto &atom : atom_list_) {
    auto cartesian_position = atom.GetCartesianPosition();
    for (const auto kDim : All_Dimensions) {
      double displacement;
      // make sure displacement is not too large
      do {
        displacement = distribution(generator);
      } while (std::abs(displacement) > kPerturbCutOff);
      cartesian_position[kDim] += displacement;
    }
    atom.SetCartesianPosition(cartesian_position);
  }
  WrapAtomCartesian();
}
void Config::ClearNeighbors() {
  for (auto &atom:atom_list_)
    atom.CleanNeighborsLists();
}
// TODO rewrite this function
void Config::UpdateNeighbors(double first_r_cutoff,
                             double second_r_cutoff,
                             double third_r_cutoff) {
  ClearNeighbors();

  const double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  const double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  const double third_r_cutoff_square = third_r_cutoff * third_r_cutoff;

  for (auto it1 = atom_list_.begin(); it1 != atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 != it1; ++it2) {
      Vector_t absolute_distance_vector = GetRelativeDistanceVector(*it1, *it2) * basis_;
      if (std::abs(absolute_distance_vector[kXDimension]) > third_r_cutoff)
        continue;
      if (std::abs(absolute_distance_vector[kYDimension]) > third_r_cutoff)
        continue;
      if (std::abs(absolute_distance_vector[kZDimension]) > third_r_cutoff)
        continue;
      double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square <= third_r_cutoff_square) {
        if (absolute_distance_square <= second_r_cutoff_square) {
          if (absolute_distance_square <= first_r_cutoff_square) {
            it1->AppendFirstNearestNeighborsList(it2->GetId());
            it2->AppendFirstNearestNeighborsList(it1->GetId());
          } else {
            it1->AppendSecondNearestNeighborsList(it2->GetId());
            it2->AppendSecondNearestNeighborsList(it1->GetId());
          }
        } else {
          it1->AppendThirdNearestNeighborsList(it2->GetId());
          it2->AppendThirdNearestNeighborsList(it1->GetId());
        }
      }
    }
  }
}
void Config::AppendAtomWithoutChangingAtomID(const Atom &atom) {
  atom_list_.push_back(atom);
  // element_list_map_[atom.GetType()].emplace_back(atom.GetId());
}
void Config::AppendAtomWithChangingAtomID(Atom atom) {
  atom.SetId(atom_list_.size());
  // element_list_map_[atom.GetType()].emplace_back(atom.GetId());
  atom_list_.push_back(std::move(atom));
}
void Config::RemoveAtomWithID(size_t id) {
  atom_list_.erase(std::find_if(atom_list_.begin(), atom_list_.end(), [id](const auto &atom) {
    return atom.GetId() == id;
  }));
}
void Config::SortAtomListWith(const std::function<bool(const cfg::Atom &,
                                                       const cfg::Atom &)> &compare) {
  //sort using mm2 group point
  std::sort(atom_list_.begin(), atom_list_.end(), compare);
}
void Config::ReIndexAtoms() {
  size_t new_id = 0;
  for (auto &atom : atom_list_) {
    atom.SetId(new_id++);
  }
}
void Config::ChangeAtomTypeAt(size_t id, const std::string &type) {
  atom_list_[id].SetType(type);
}

Config Config::ReadPOSCAR(const std::string &filename, bool update_neighbors) {
  std::ifstream ifs(filename, std::ifstream::in);
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // #comment
  double scale;
  ifs >> scale; // scale factor, usually which is 1.0
  Matrix_t basis;
  ifs >> basis;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line

  std::string buffer;
  getline(ifs, buffer);
  std::istringstream element_iss(buffer);
  getline(ifs, buffer);
  std::istringstream count_iss(buffer);

  std::string element;
  size_t num_atoms;
  size_t all_num_atoms = 0;
  std::vector<std::pair<std::string, size_t>> elements_counts;
  while (element_iss >> element && count_iss >> num_atoms) {
    elements_counts.emplace_back(element, num_atoms);
    all_num_atoms += num_atoms;
  }
  getline(ifs, buffer);
  bool relative_option =
      buffer[0] != 'C' && buffer[0] != 'c' && buffer[0] != 'K' && buffer[0] != 'k';
  Config config(basis * scale, all_num_atoms);

  // If relative_option is ture, only relative position need to scaled, set it to 1
  if (relative_option)
    scale = 1.0;

  size_t id_count = 0;
  double position_X, position_Y, position_Z;
  for (const auto &[element_name, count] : elements_counts) {
    double mass = FindMass(element_name);
    for (size_t j = 0; j < count; ++j) {
      ifs >> position_X >> position_Y >> position_Z;
      config.AppendAtomWithoutChangingAtomID({id_count, mass, element_name,
                                              position_X * scale, position_Y * scale,
                                              position_Z * scale});
      ++id_count;
    }
  }
  if (relative_option)
    config.ConvertRelativeToCartesian();
  else
    config.ConvertCartesianToRelative();

  if (update_neighbors)
    config.UpdateNeighbors();
  return config;
}
Config Config::ReadConfig(const std::string &filename, bool update_neighbors) {
  std::ifstream ifs(filename, std::ifstream::in);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "Number of particles = %i"
  size_t num_atoms;
  ifs >> num_atoms;

  ifs.ignore(std::numeric_limits<std::streamsize>::max(),
             '='); // A = 1.0 Angstrom (basic length-scale)
  double scale;
  ifs >> scale;

  double basis_xx, basis_xy, basis_xz, basis_yx, basis_yy, basis_yz, basis_zx, basis_zy,
      basis_zz;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(1,1) = %lf A"
  ifs >> basis_xx;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(1,2) = %lf A"
  ifs >> basis_xy;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(1,3) = %lf A"
  ifs >> basis_xz;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(2,1) = %lf A"
  ifs >> basis_yx;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(2,2) = %lf A"
  ifs >> basis_yy;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(2,3) = %lf A"
  ifs >> basis_yz;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(3,1) = %lf A"
  ifs >> basis_zx;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(3,2) = %lf A"
  ifs >> basis_zy;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "H0(3,3) = %lf A"
  ifs >> basis_zz;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line
  Config config(Matrix_t{{{basis_xx, basis_xy, basis_xz},
                          {basis_yx, basis_yy, basis_yz},
                          {basis_zx, basis_zy, basis_zz}}} * scale, num_atoms);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // .NO_VELOCITY.
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // "entry_count = 3"

  double mass, relative_position_X, relative_position_Y, relative_position_Z;
  size_t index;
  bool neighbor_found = false;
  for (size_t id = 0; id < num_atoms; ++id) {
    std::string type;
    ifs >> mass;
    ifs >> type;
    ifs >> relative_position_X >> relative_position_Y >> relative_position_Z;
    Atom atom(id, mass, type,
              relative_position_X * scale,
              relative_position_Y * scale,
              relative_position_Z * scale);
    if (ifs.peek() != '\n') {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '#');
      for (size_t i = 0; i < Al_const::kNumFirstNearestNeighbors; ++i) {
        ifs >> index;
        atom.AppendFirstNearestNeighborsList(index);
      }
      for (size_t i = 0; i < Al_const::kNumSecondNearestNeighbors; ++i) {
        ifs >> index;
        atom.AppendSecondNearestNeighborsList(index);
      }
      for (size_t i = 0; i < Al_const::kNumThirdNearestNeighbors; ++i) {
        ifs >> index;
        atom.AppendThirdNearestNeighborsList(index);
      }
      neighbor_found = true;
    }
    config.AppendAtomWithoutChangingAtomID(atom);
  }
  config.ConvertRelativeToCartesian();
  if (!neighbor_found && update_neighbors)
    config.UpdateNeighbors();
  return config;
}
// Write Configuration out as POSCAR file. If the show_vacancy_option is
// true, output will have "X" for visualization. If false, vacancies will be
// ignored for VASP calculation.
void Config::WritePOSCAR(const Config &config,
                         const std::string &filename,
                         bool show_vacancy_option) {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << "#comment\n1.0\n";
  ofs << config.GetBasis() << '\n';
  std::ostringstream ele_oss, count_oss;
  const auto element_list_map = config.GetElementListMap();
  for (const auto &[element, element_list] : element_list_map) {
    if (!show_vacancy_option && element == "X")
      continue;
    ele_oss << element << ' ';
    count_oss << element_list.size() << ' ';
  }
  ofs << ele_oss.str() << '\n' << count_oss.str() << '\n';
  ofs << "Direct\n";

  for (const auto &[element, element_list] : element_list_map) {
    if (!show_vacancy_option && element == "X")
      continue;
    for (auto index : element_list) {
      ofs << config.GetAtomList()[index].GetRelativePosition() << '\n';
    }
  }
}
void Config::WriteConfig(const Config &config,
                         const std::string &filename,
                         bool neighbors_info) {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs.precision(16);
  ofs << "Number of particles = " << config.GetNumAtoms() << '\n';
  ofs << "A = 1.0 Angstrom (basic length-scale)\n";
  auto basis = config.GetBasis();
  ofs << "H0(1,1) = " << basis[kXDimension][kXDimension] << " A\n";
  ofs << "H0(1,2) = " << basis[kXDimension][kYDimension] << " A\n";
  ofs << "H0(1,3) = " << basis[kXDimension][kZDimension] << " A\n";
  ofs << "H0(2,1) = " << basis[kYDimension][kXDimension] << " A\n";
  ofs << "H0(2,2) = " << basis[kYDimension][kYDimension] << " A\n";
  ofs << "H0(2,3) = " << basis[kYDimension][kZDimension] << " A\n";
  ofs << "H0(3,1) = " << basis[kZDimension][kXDimension] << " A\n";
  ofs << "H0(3,2) = " << basis[kZDimension][kYDimension] << " A\n";
  ofs << "H0(3,3) = " << basis[kZDimension][kZDimension] << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom : config.GetAtomList()) {
    ofs << atom.GetMass() << '\n'
        << atom.GetType() << '\n'
        << atom.GetRelativePosition();
    if (neighbors_info) {
      ofs << " # ";
      for (auto neighbor_index : atom.GetFirstNearestNeighborsList()) {
        ofs << neighbor_index << ' ';
      }
      for (auto neighbor_index : atom.GetSecondNearestNeighborsList()) {
        ofs << neighbor_index << ' ';
      }
      for (auto neighbor_index : atom.GetThirdNearestNeighborsList()) {
        ofs << neighbor_index << ' ';
      }
    }
    ofs << '\n';
    ofs << std::flush;
  }
}
std::unordered_set<size_t> GetFirstAndSecondThirdNeighborsSetOfJumpPair(
    const Config &config, const std::pair<size_t, size_t> &jump_pair) {

  std::unordered_set<size_t> near_neighbors_hashset;
  for (const auto i : {jump_pair.first, jump_pair.second}) {
    const auto &atom = config.GetAtomList()[i];
    std::copy(atom.GetFirstNearestNeighborsList().begin(),
              atom.GetFirstNearestNeighborsList().end(),
              std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
    std::copy(atom.GetSecondNearestNeighborsList().begin(),
              atom.GetSecondNearestNeighborsList().end(),
              std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
    std::copy(atom.GetThirdNearestNeighborsList().begin(),
              atom.GetThirdNearestNeighborsList().end(),
              std::inserter(near_neighbors_hashset, near_neighbors_hashset.begin()));
  }
  return near_neighbors_hashset;
}

// lhs : vacancy_index; rhs : neighbor_index
void AtomsJump(Config &config, const std::pair<size_t, size_t> &jump_pair) {
  const auto[lhs, rhs] = jump_pair;
  auto &atom_list = config.atom_list_;
  //   auto it_lhs = std::find_if(config.atom_list_.begin(),
  //                              config.atom_list_.end(),
  //                              [jump_pair](const auto &atom) {
  //                                return atom.GetId() == jump_pair.first;
  //                              });
  //   auto it_rhs = std::find_if(config.atom_list_.begin(),
  //                              config.atom_list_.end(),
  //                              [jump_pair](const auto &atom) {
  //                                return atom.GetId() == jump_pair.second;
  //                              });
  //   const auto lhs = static_cast<const size_t>(std::distance(config.atom_list_.begin(), it_lhs));
  //   const auto rhs = static_cast<const size_t>(std::distance(config.atom_list_.begin(), it_rhs));

  // three things to swap here:
  // 1) element coordinates
  std::swap(atom_list[lhs].relative_position_,
            atom_list[rhs].relative_position_);
  std::swap(atom_list[lhs].cartesian_position_,
            atom_list[rhs].cartesian_position_);

  // 2) jump pair neighbour lists
  std::swap(atom_list[lhs].first_nearest_neighbors_list_,
            atom_list[rhs].first_nearest_neighbors_list_);

  std::swap(atom_list[lhs].second_nearest_neighbors_list_,
            atom_list[rhs].second_nearest_neighbors_list_);

  std::swap(atom_list[lhs].third_nearest_neighbors_list_,
            atom_list[rhs].third_nearest_neighbors_list_);

  // 3) jump atoms' and neighbor atom's neighbor list
  std::unordered_set<size_t>
      atom_id_hashset = GetFirstAndSecondThirdNeighborsSetOfJumpPair(config, jump_pair);

  for (auto i : atom_id_hashset) {
    // auto it_id = std::find_if(config.atom_list_.begin(),
    //                           config.atom_list_.end(),
    //                           [id](const auto &atom) {
    //                             return atom.GetId() == id;
    //                           });
    // const auto i = static_cast<const size_t>(std::distance(config.atom_list_.begin(), it_id));
    for (auto neighbors_list : {&atom_list[i].first_nearest_neighbors_list_,
                                &atom_list[i].second_nearest_neighbors_list_,
                                &atom_list[i].third_nearest_neighbors_list_})
      for (auto &j : *neighbors_list) {
        if (j == lhs)
          j = rhs;
        else if (j == rhs)
          j = lhs;
      }
  }
}
void AtomsJumpMore(Config &config, const std::pair<size_t, size_t> &jump_pair) {
  const auto[lhs, rhs] = jump_pair;
  auto &atom_list = config.atom_list_;
  // three things to swap here:
  // 1) element coordinates
  std::swap(atom_list[lhs].relative_position_,
            atom_list[rhs].relative_position_);
  std::swap(atom_list[lhs].cartesian_position_,
            atom_list[rhs].cartesian_position_);

  // 2) jump pair neighbour lists
  std::swap(atom_list[lhs].first_nearest_neighbors_list_,
            atom_list[rhs].first_nearest_neighbors_list_);

  std::swap(atom_list[lhs].second_nearest_neighbors_list_,
            atom_list[rhs].second_nearest_neighbors_list_);

  std::swap(atom_list[lhs].third_nearest_neighbors_list_,
            atom_list[rhs].third_nearest_neighbors_list_);

  // 3) jump atoms' and neighbor atom's neighbor list
  std::unordered_set<size_t>
      atom_id_hashset = GetFirstAndSecondThirdNeighborsSetOfJumpPair(config, jump_pair);

  for (auto i : atom_id_hashset) {
    for (auto neighbors_list : {&atom_list[i].first_nearest_neighbors_list_,
                                &atom_list[i].second_nearest_neighbors_list_,
                                &atom_list[i].third_nearest_neighbors_list_})
      for (auto &j : *neighbors_list) {
        if (j == lhs)
          j = rhs;
        else if (j == rhs)
          j = lhs;
      }
  }

}
std::map<std::string, size_t> CountAllType(const Config &config) {
  std::map<std::string, size_t> types_count_map;
  for (const auto &atom : config.GetAtomList()) {
    types_count_map[atom.GetType()]++;
  }
  return types_count_map;
}

std::map<Bond, size_t> CountAllBonds(const Config &config) {
  std::map<Bond, size_t> bonds_count_map;
  std::string type1, type2;
  const auto &atom_list = config.GetAtomList();
  for (const auto &atom : atom_list) {
    type1 = atom.GetType();
    for (const auto &atom2_id : atom.GetFirstNearestNeighborsList()) {
      bonds_count_map[Bond{type1, atom_list[atom2_id].GetType()}]++;
    }
  }
  for (auto &bond_count : bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
std::unordered_map<std::string, size_t> GetTypeCategoryHashmap(const Config &config) {
  size_t count = 1;
  std::unordered_map<std::string, size_t> type_category_hashmap;
  for (const auto &element_list : config.GetElementListMap()) {
    type_category_hashmap[element_list.first] = count++;
  }
  return type_category_hashmap;
}

Vector_t GetPairCenter(const Config &config, const std::pair<size_t, size_t> &jump_pair) {
  Vector_t center_position;
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
Matrix_t GetPairRotationMatrix(const Config &config, const std::pair<size_t, size_t> &jump_pair) {
  const Vector_t pair_direction = Normalize(GetRelativeDistanceVector(
      config.GetAtomList()[jump_pair.first],
      config.GetAtomList()[jump_pair.second]));
  const auto &first_atom = config.GetAtomList()[jump_pair.first];
  Vector_t vertical_vector{};
  for (const auto index : first_atom.GetFirstNearestNeighborsList()) {
    const Vector_t jump_vector = GetRelativeDistanceVector(first_atom, config.GetAtomList()[index]);
    const double dot_prod = Dot(pair_direction, jump_vector);
    if (std::abs(dot_prod) < 1e-6) {
      vertical_vector = Normalize(jump_vector);
      break;
    }
  }
  // The third row is normalized since it is a cross product of two normalized vectors.
  // We use transposed matrix here because transpose of an orthogonal matrix equals its inverse
  return TransposeMatrix({pair_direction, vertical_vector,
                          Cross(pair_direction, vertical_vector)});
}
void RotateAtomVector(std::vector<Atom> &atom_list, const Matrix_t &rotation_matrix) {
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
}
size_t GetVacancyIndex(const Config &config) {
  const auto &atom_vector = config.GetAtomList();
  auto it = std::find_if(atom_vector.cbegin(), atom_vector.cend(), [](const auto &atom) {
    return atom.GetType() == "X";
  });
  if (it != atom_vector.end()) {
    return it->GetId();
  } else {
    std::cout << "Warning: vacancy not found" << std::endl;
  }
  return 0;
}
Config GenerateFCC(double lattice_constant_a, const std::string &element, const Factor_t &factors) {

  double mass = FindMass(element);
  Config config({{{lattice_constant_a * static_cast<double>(factors[kXDimension]), 0, 0},
                  {0, lattice_constant_a * static_cast<double>(factors[kYDimension]), 0},
                  {0, 0, lattice_constant_a * static_cast<double>(factors[kZDimension])}}},
                4 * factors[kXDimension] * factors[kYDimension] * factors[kZDimension]);
  size_t atoms_counter = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (size_t k = 0; k < factors[kZDimension]; ++k) {
    for (size_t j = 0; j < factors[kYDimension]; ++j) {
      for (size_t i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        config.AppendAtomWithoutChangingAtomID({
                                                   atoms_counter++, mass, element,
                                                   x_reference / x_length,
                                                   y_reference / y_length,
                                                   z_reference / z_length
                                               });
        config.AppendAtomWithoutChangingAtomID({
                                                   atoms_counter++, mass, element,
                                                   (x_reference + 0.5) / x_length,
                                                   (y_reference + 0.5) / y_length,
                                                   z_reference / z_length
                                               });
        config.AppendAtomWithoutChangingAtomID({
                                                   atoms_counter++, mass, element,
                                                   (x_reference + 0.5) / x_length,
                                                   y_reference / y_length,
                                                   (z_reference + 0.5) / z_length
                                               });
        config.AppendAtomWithoutChangingAtomID({
                                                   atoms_counter++, mass, element,
                                                   x_reference / x_length,
                                                   (y_reference + 0.5) / y_length,
                                                   (z_reference + 0.5) / z_length
                                               });
      }
    }
  }
  config.ConvertRelativeToCartesian();
  config.UpdateNeighbors();
  return config;
}

}// namespace cfg