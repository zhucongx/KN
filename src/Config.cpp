#include"Config.h"

#include <iostream>
#include <sstream>
#include <utility>
#include <random>
#include <chrono>

namespace box {

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kPerturbCutOff = 0.4;

Config::Config() = default;

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}

void Config::Initialize() {
  basis_ = {{{0, 0, 0},
             {0, 0, 0},
             {0, 0, 0}}};
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_map_.clear();
  neighbor_found_ = false;
}

bool Config::IsCubic() const {
  return basis_[kXDimension][kXDimension] == basis_[kYDimension][kYDimension] &&
      basis_[kYDimension][kYDimension] == basis_[kZDimension][kZDimension] &&
      basis_[kXDimension][kYDimension] == 0 &&
      basis_[kXDimension][kZDimension] == 0 &&
      basis_[kYDimension][kXDimension] == 0 &&
      basis_[kYDimension][kZDimension] == 0 &&
      basis_[kZDimension][kXDimension] == 0 &&
      basis_[kZDimension][kYDimension] == 0;
}

void Config::ConvertRelativeToCartesian() {
  for (auto &atom:atom_list_) {
    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}

void Config::ConvertCartesianToRelative() {
  auto inverse_basis = InverseMatrix33(basis_);
  for (auto &atom:atom_list_) {
    atom.relative_position_ = atom.cartesian_position_ * inverse_basis;
  }
}

void Config::Perturb() {
  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937_64 generator(seed);
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  auto add_displacement = [&generator, &distribution](double &coordinate) {
    double displacement = distribution(generator);
    while (std::abs(displacement) > kPerturbCutOff) {
      displacement = distribution(generator);
    }
    coordinate += displacement;
  };
  for (auto &atom:atom_list_) {
    add_displacement(atom.cartesian_position_[kXDimension]);
    add_displacement(atom.cartesian_position_[kYDimension]);
    add_displacement(atom.cartesian_position_[kZDimension]);
  }
  ConvertCartesianToRelative();
}

void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff) {
  if (neighbor_found_)
    return;

  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Vector3 absolute_distance_vector =
          GetRelativeDistanceVector(*it1, *it2) * basis_;
      if (absolute_distance_vector[kXDimension] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kYDimension] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kZDimension] > second_r_cutoff_square)
        continue;
      double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square <= second_r_cutoff_square) {
        if (absolute_distance_square <= first_r_cutoff_square) {
          it1->first_nearest_neighbor_list_.emplace_back(it2->GetId());
          it2->first_nearest_neighbor_list_.emplace_back(it1->GetId());
        }
        it1->second_nearest_neighbor_list_.emplace_back(it2->GetId());
        it2->second_nearest_neighbor_list_.emplace_back(it1->GetId());
      }
    }
  }
  neighbor_found_ = true;
}

void Config::WrapRelativePosition() {
  for (auto &atom:atom_list_) {
    atom.relative_position_[kXDimension] -= floor(atom.relative_position_[kXDimension]);
    atom.relative_position_[kYDimension] -= floor(atom.relative_position_[kYDimension]);
    atom.relative_position_[kZDimension] -= floor(atom.relative_position_[kZDimension]);

    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}

// void Config::WrapAbsolutePosition() {
//   ConvertCartesianToRelative();
//   WrapRelativePosition();
// }

// void Config::ShiftAtomToCentral(const Atom::Rank &id) {
//   Vector3 criterionDistance =
//       atom_list_[id].relative_position_ - Vector3{0.5, 0.5, 0.5};
//   for (auto &atom:atom_list_) {
//     atom.relative_position_ = atom.relative_position_ - criterionDistance;
//   }
//   WrapRelativePosition();
// }
// for better performance, shouldn't call Wrap function
void Config::MoveRelativeDistance(const Vector3 &distance_vector) {
  for (auto &atom:atom_list_) {
    atom.relative_position_ += distance_vector;

    atom.relative_position_ -= ElementFloor(atom.relative_position_);

    atom.cartesian_position_ = atom.relative_position_ * basis_;
  }
}
void Config::MoveOneAtomRelativeDistance(const Atom::Rank &index,
                                         const Vector3 &distance_vector) {
  atom_list_[index].relative_position_ += distance_vector;
  atom_list_[index].relative_position_ -= ElementFloor(atom_list_[index].relative_position_);

  atom_list_[index].cartesian_position_ = atom_list_[index].relative_position_ * basis_;
}
// void Config::MoveAbsoluteDistance(const Vector3 &distance_vector) {
//   for (auto &atom:atom_list_) {
//     atom.cartesian_position_ = atom.cartesian_position_ + distance_vector;
//   }
//   WrapAbsolutePosition();
// }

std::map<Bond, int> Config::CountAllBonds(double r_cutoff) {
  if (!neighbor_found_)
    UpdateNeighbors(r_cutoff, r_cutoff);

  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  for (const auto &atom:atom_list_) {
    type1 = atom.GetType();
    for (const auto &atom2_id:atom.first_nearest_neighbor_list_) {
      type2 = atom_list_[atom2_id].GetType();
      bonds_count_map[Bond{type1, type2}]++;
    }
  }
  for (auto &bond_count:bonds_count_map) {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
void Config::ReadConfig(const std::string &file_name) {
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "Number of particles = %i"
  int num_atoms;
  ifs >> num_atoms;

  ifs.ignore(std::numeric_limits<std::streamsize>::max(),
             '=');  // A = 1.0 Angstrom (basic length-scale)
  ifs >> scale_;

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(1,1) = %lf A"
  ifs >> basis_[kXDimension][kXDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(1,2) = %lf A"
  ifs >> basis_[kXDimension][kYDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(1,3) = %lf A"
  ifs >> basis_[kXDimension][kZDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(2,1) = %lf A"
  ifs >> basis_[kYDimension][kXDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(2,2) = %lf A"
  ifs >> basis_[kYDimension][kYDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(2,3) = %lf A"
  ifs >> basis_[kYDimension][kZDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(3,1) = %lf A"
  ifs >> basis_[kZDimension][kXDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(3,2) = %lf A"
  ifs >> basis_[kZDimension][kYDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '=');  // "H0(3,3) = %lf A"
  ifs >> basis_[kZDimension][kZDimension];

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // .NO_VELOCITY.
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // "entry_count = 3"
  for (Atom::Rank i = 0; i < num_atoms; ++i) {
    double mass, relative_position_X, relative_position_Y, relative_position_Z;
    std::string type;
    ifs >> mass;
    ifs >> type;
    ifs >> relative_position_X >> relative_position_Y >> relative_position_Z;
    atom_list_.emplace_back(i, mass, type,
                            relative_position_X,
                            relative_position_Y,
                            relative_position_Z);
    element_list_map_[type].emplace_back(i);
  }
  ConvertRelativeToCartesian();
}

void Config::ReadPOSCAR(const std::string &file_name) {
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // #comment

  ifs >> scale_;  // scale factor, usually which is 1.0
  ifs >> basis_;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // finish this line

  std::string buffer;
  getline(ifs, buffer);
  std::istringstream element_iss(buffer);
  getline(ifs, buffer);
  std::istringstream count_iss(buffer);

  std::string element;
  int count;

  std::vector<std::pair<std::string, int>> elements_counts;
  while (element_iss >> element && count_iss >> count) {
    elements_counts.emplace_back(element, count);
  }
  getline(ifs, buffer);
  bool relative_option;
  relative_option = buffer[0] == 'D' || buffer[0] == 'd';

  Atom::Rank id_count = 0;
  double position_X, position_Y, position_Z;
  for (const auto&[element_name, count]:elements_counts) {
    double mass = elem_info::FindMass(element_name);
    for (int j = 0; j < count; ++j) {
      ifs >> position_X >> position_Y >> position_Z;
      atom_list_.emplace_back(id_count, mass, element_name,
                              position_X, position_Y, position_Z);
      element_list_map_[element_name].emplace_back(id_count);
      ++id_count;
    }
  }
  if (relative_option) {
    ConvertRelativeToCartesian();
  } else {
    ConvertCartesianToRelative();
  }
}

void Config::WriteConfig(const std::string &file_name) const {
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "Number of particles = " << atom_list_.size() << '\n';
  ofs << "A = " << scale_ << " Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << basis_[kXDimension][kXDimension] << " A\n";
  ofs << "H0(1,2) = " << basis_[kXDimension][kYDimension] << " A\n";
  ofs << "H0(1,3) = " << basis_[kXDimension][kZDimension] << " A\n";
  ofs << "H0(2,1) = " << basis_[kYDimension][kXDimension] << " A\n";
  ofs << "H0(2,2) = " << basis_[kYDimension][kYDimension] << " A\n";
  ofs << "H0(2,3) = " << basis_[kYDimension][kZDimension] << " A\n";
  ofs << "H0(3,1) = " << basis_[kZDimension][kXDimension] << " A\n";
  ofs << "H0(3,2) = " << basis_[kZDimension][kYDimension] << " A\n";
  ofs << "H0(3,3) = " << basis_[kZDimension][kZDimension] << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom:atom_list_) {
    double mass = atom.GetMass();
    const std::string &type = atom.GetType();
    ofs << mass << '\n'
        << type << '\n'
        << atom.relative_position_ << '\n';
  }
}

void Config::WritePOSCAR(const std::string &file_name,
                         const bool &show_vacancy_option) const {
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "#comment\n" << scale_ << '\n';
  ofs << basis_ << '\n';
  std::ostringstream ele_oss, count_oss;
  for (const auto &[element, element_list]:element_list_map_) {
    if (!show_vacancy_option || element != "X") {
      ele_oss << element << ' ';
      count_oss << element_list.size() << ' ';
    }
  }
  ofs << ele_oss.str() << '\n' << count_oss.str() << '\n';
  ofs << "Direct\n";
  for (const auto &atom:atom_list_) {
    if (!show_vacancy_option || atom.GetType() != "X") {
      ofs << atom.relative_position_ << '\n';
    }
  }
}
void Config::GenerateUnitCell(const Matrix33 &bravais_matrix,
                              const std::vector<std::pair<std::string,
                                                          Vector3>> &type_position_list) {
  Initialize();
  basis_ = bravais_matrix;
  scale_ = 1.0;
  int atoms_counter_ = 0;
  for (const auto&[type, relative_position]:type_position_list) {
    atom_list_.emplace_back(atoms_counter_,
                            elem_info::FindMass(type),
                            type,
                            relative_position);
    element_list_map_[type].emplace_back(atoms_counter_++);
  }
  ConvertRelativeToCartesian();
}
void Config::Duplicate(const std::array<int, kDimension> &factors) {
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  basis_ = {basis_[kXDimension] * x_length,
            basis_[kYDimension] * y_length,
            basis_[kZDimension] * z_length};
  auto temp(std::move(atom_list_));
  int atoms_counter_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_map_.clear();
  neighbor_found_ = false;
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        for (const auto &atom:temp) {
          atom_list_.emplace_back(atoms_counter_, atom.GetMass(), atom.GetType(),
                                  (x_reference + atom.relative_position_[kXDimension])
                                      / x_length,
                                  (y_reference + atom.relative_position_[kYDimension])
                                      / y_length,
                                  (z_reference + atom.relative_position_[kZDimension])
                                      / z_length);
          element_list_map_[atom.GetType()].emplace_back(atoms_counter_++);
        }
      }
    }
  }

}
void Config::GenerateFCC(const double &lattice_constant_a,
                         const std::string &element,
                         const std::array<int, kDimension> &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  basis_ = {{{lattice_constant_a * factors[kXDimension], 0, 0},
             {0, lattice_constant_a * factors[kYDimension], 0},
             {0, 0, lattice_constant_a * factors[kZDimension]}}};
  scale_ = 1.0;
  int atoms_counter_ = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                z_reference / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                (x_reference + 0.5) / x_length,
                                y_reference / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                x_reference / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
      }
    }
  }
  ConvertRelativeToCartesian();
}
void Config::GenerateBCC(const double &lattice_constant_a,
                         const std::string &element,
                         const std::array<int, kDimension> &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  basis_ = {{{lattice_constant_a * factors[kXDimension], 0, 0},
             {0, lattice_constant_a * factors[kYDimension], 0},
             {0, 0, lattice_constant_a * factors[kZDimension]}}};
  scale_ = 1.0;
  int atoms_counter_ = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
      }
    }
  }
  ConvertRelativeToCartesian();
}
void Config::GenerateHCP(const double &lattice_constant_a,
                         const double &lattice_constant_c,
                         const std::string &element,
                         const std::array<int, kDimension> &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  basis_ = {{{lattice_constant_a * factors[kXDimension], 0, 0},
             {-0.5 * lattice_constant_a * factors[kYDimension],
              0.5 * sqrt(3) * lattice_constant_a * factors[kYDimension], 0},
             {0, 0, lattice_constant_c * factors[kZDimension]}}};
  scale_ = 1.0;
  int atoms_counter_ = 0;

  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k) {
    for (int j = 0; j < factors[kYDimension]; ++j) {
      for (int i = 0; i < factors[kXDimension]; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
        atom_list_.emplace_back(atoms_counter_, mass, element,
                                (x_reference + 1.0 / 3.0) / x_length,
                                (y_reference + 2.0 / 3.0) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_map_[element].emplace_back(atoms_counter_++);
      }
    }
  }
  ConvertRelativeToCartesian();
}
void Config::SetScale(double scale) {
  scale_ = scale;
}

double Config::GetScale() const {
  return scale_;
}
const Matrix33 &Config::GetBasis() const {
  return basis_;
}
void Config::SetBasis(const Matrix33 &basis) {
  basis_ = basis;
}
const Atom &Config::GetAtom(const int &index) const {
  return atom_list_[index];
}
void Config::AppendAtom(const Atom &atom) {
  atom_list_.push_back(atom);
  element_list_map_[atom.GetType()].emplace_back(atom.GetId());
}
int Config::GetNumAtoms() const {
  return atom_list_.size();
}
const std::map<std::string, std::vector<Atom::Rank>> &Config::GetElementListMap() const {
  return element_list_map_;
}

}// namespace box
