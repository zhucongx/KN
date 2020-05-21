#include"Config.h"

#include <iostream>
#include <sstream>
#include <utility>
#include <random>
#include <chrono>

namespace box
{

const double kMean = 0;
const double kStandardDeviation = 0.15;
const double kPerturbCutOff = 0.4;

Config::Config() = default;

bool Config::operator<(const Config &rhs) const
{
  return energy_ < rhs.energy_;
}

void Config::Initialize()
{
  bravais_matrix_ = {{{0, 0, 0},
                      {0, 0, 0},
                      {0, 0, 0}}};
  num_atoms_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_set_.clear();
  neighbor_found_ = false;
}

bool Config::IsCubic() const
{
  return bravais_matrix_[kXDimension][kXDimension] == bravais_matrix_[kYDimension][kYDimension] &&
      bravais_matrix_[kYDimension][kYDimension] == bravais_matrix_[kZDimension][kZDimension] &&
      bravais_matrix_[kXDimension][kYDimension] == 0 &&
      bravais_matrix_[kXDimension][kZDimension] == 0 &&
      bravais_matrix_[kYDimension][kXDimension] == 0 &&
      bravais_matrix_[kYDimension][kZDimension] == 0 &&
      bravais_matrix_[kZDimension][kXDimension] == 0 &&
      bravais_matrix_[kZDimension][kYDimension] == 0;
}

void Config::ConvertRelativeToAbsolute()
{
  for (auto &atom:atom_list_)
  {
    atom.cartesian_position_ = atom.relative_position_ * bravais_matrix_;
  }
}

void Config::ConvertAbsoluteToRelative()
{
  for (auto &atom:atom_list_)
  {
    atom.relative_position_ = atom.cartesian_position_ * inverse_bravais_matrix_;
  }
}

void Config::Perturb()
{
  auto seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937_64 generator(seed);
  std::normal_distribution<double> distribution(kMean, kStandardDeviation);
  auto add_displacement = [&generator, &distribution](double &coordinate)
  {
    double displacement = distribution(generator);
    while (std::abs(displacement) > kPerturbCutOff)
    {
      displacement = distribution(generator);
    }
    coordinate += displacement;
  };
  for (auto &atom:atom_list_)
  {
    add_displacement(atom.cartesian_position_[kXDimension]);
    add_displacement(atom.cartesian_position_[kYDimension]);
    add_displacement(atom.cartesian_position_[kZDimension]);
  }
  ConvertAbsoluteToRelative();
}

void Config::UpdateNeighbors(double first_r_cutoff, double second_r_cutoff)
{
  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;
  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1)
  {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2)
    {
      Vector3 absolute_distance_vector =
          GetRelativeDistanceVector(*it1, *it2) * bravais_matrix_;
      if (absolute_distance_vector[kXDimension] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kYDimension] > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector[kZDimension] > second_r_cutoff_square)
        continue;
      double absolute_distance_square = Inner(absolute_distance_vector);
      if (absolute_distance_square <= second_r_cutoff_square)
      {
        if (absolute_distance_square <= first_r_cutoff_square)
        {
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

void Config::WrapRelativePosition()
{
  for (auto &atom:atom_list_)
  {
    atom.relative_position_[kXDimension] -= floor(atom.relative_position_[kXDimension]);
    atom.relative_position_[kYDimension] -= floor(atom.relative_position_[kYDimension]);
    atom.relative_position_[kZDimension] -= floor(atom.relative_position_[kZDimension]);

    atom.cartesian_position_ = atom.relative_position_ * bravais_matrix_;
  }
}

// void Config::WrapAbsolutePosition() {
//   ConvertAbsoluteToRelative();
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
void Config::MoveRelativeDistance(const Vector3 &distance_vector)
{
  for (auto &atom:atom_list_)
  {
    atom.relative_position_ += distance_vector;

    atom.relative_position_ -= ElementFloor(atom.relative_position_);

    atom.cartesian_position_ = atom.relative_position_ * bravais_matrix_;
  }
}
void Config::MoveOneAtomRelativeDistance(const Atom::Rank &index,
                                         const Vector3 &distance_vector)
{
  atom_list_[index].relative_position_ += distance_vector;
  atom_list_[index].relative_position_ -= ElementFloor(atom_list_[index].relative_position_);

  atom_list_[index].cartesian_position_ = atom_list_[index].relative_position_ * bravais_matrix_;
}
// void Config::MoveAbsoluteDistance(const Vector3 &distance_vector) {
//   for (auto &atom:atom_list_) {
//     atom.cartesian_position_ = atom.cartesian_position_ + distance_vector;
//   }
//   WrapAbsolutePosition();
// }

std::map<Bond, int> Config::CountAllBonds(double r_cutoff)
{
  if (!neighbor_found_)
    UpdateNeighbors(r_cutoff, r_cutoff);

  std::map<Bond, int> bonds_count_map;
  std::string type1, type2;
  for (const auto &atom:atom_list_)
  {
    type1 = atom.GetType();
    for (const auto &atom2_id:atom.first_nearest_neighbor_list_)
    {
      type2 = atom_list_[atom2_id].GetType();
      bonds_count_map[Bond{type1, type2}]++;
    }
  }
  for (auto &bond_count:bonds_count_map)
  {
    bond_count.second /= 2;
  }
  return bonds_count_map;
}
bool Config::ReadConfig(const std::string &file_name)
{
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);
  if (ifs.fail())
  {
    return false;
  }
  std::string line;
  std::istringstream iss;
  if (!getline(ifs, line))  // "Number of particles = %i"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> num_atoms_))
  {
    return false;
  }
  if (!getline(ifs, line))  // A = 1.0 Angstrom (basic length-scale)
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> scale_))
  {
    return false;
  }

  if (!getline(ifs, line))  // "H0(1,1) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kXDimension][kXDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(1,2) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kXDimension][kYDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(1,3) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kXDimension][kZDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(2,1) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kYDimension][kXDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(2,2) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kYDimension][kYDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(2,3) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kYDimension][kZDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(3,1) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kZDimension][kXDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(3,2) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kZDimension][kYDimension]))
  {
    return false;
  }
  if (!getline(ifs, line))  // "H0(3,3) = %lf A"
  {
    return false;
  }
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_[kZDimension][kZDimension]))
  {
    return false;
  }
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);

  if (!getline(ifs, line))  // .NO_VELOCITY.
  {
    return false;
  }
  if (!getline(ifs, line))  // "entry_count = 3"
  {
    return false;
  }
  for (Atom::Rank i = 0; i < num_atoms_; ++i)
  {
    double mass, relative_position_X, relative_position_Y, relative_position_Z;
    std::string type;
    if (!getline(ifs, line))
    {
      return false;
    }
    iss = std::istringstream(line);
    if (!(iss >> mass))
    {
      return false;
    }
    if (!getline(ifs, line))
    {
      return false;
    }
    type = line;
    if (!getline(ifs, line))
    {
      return false;
    }
    iss = std::istringstream(line);
    if (!(iss >> relative_position_X >> relative_position_Y
              >> relative_position_Z))
    {
      return false;
    }
    atom_list_.emplace_back(i, mass, type,
                            relative_position_X,
                            relative_position_Y,
                            relative_position_Z);
    element_list_set_[type].emplace_back(i);
  }
  ConvertRelativeToAbsolute();
  return true;
}

bool Config::ReadPOSCAR(const std::string &file_name)
{
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);
  if (ifs.fail())
  {
    return false;
  }
  std::string line;
  std::istringstream iss;
  if (!getline(ifs, line))  // #comment
  {
    return false;
  }
  if (!getline(ifs, line))  // scale factor, usually which is 1.0
  {
    return false;
  }
  iss = std::istringstream(line);
  if (!(iss >> scale_))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  iss = std::istringstream(line);
  if (!(iss >> bravais_matrix_[kXDimension][kXDimension]
            >> bravais_matrix_[kXDimension][kYDimension]
            >> bravais_matrix_[kXDimension][kZDimension]))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  iss = std::istringstream(line);
  if (!(iss >> bravais_matrix_[kYDimension][kXDimension]
            >> bravais_matrix_[kYDimension][kYDimension]
            >> bravais_matrix_[kYDimension][kZDimension]))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  iss = std::istringstream(line);
  if (!(iss >> bravais_matrix_[kZDimension][kXDimension]
            >> bravais_matrix_[kZDimension][kYDimension]
            >> bravais_matrix_[kZDimension][kZDimension]))
  {
    return false;
  }
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);

  if (!getline(ifs, line))
  {
    return false;
  }
  std::string element;
  std::istringstream name_iss(line);
  if (!getline(ifs, line))
  {
    return false;
  }
  int count;
  std::istringstream count_iss(line);

  std::vector<std::pair<std::string, int>> elements_counts;
  while (name_iss >> element && count_iss >> count)
  {
    num_atoms_ += count;
    elements_counts.emplace_back(element, count);
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  bool relOpt;
  if (line[0] == 'D' || line[0] == 'd')
  {
    relOpt = true;
  } else if (line[0] == 'C' || line[0] == 'c')
  {
    relOpt = false;
  } else
  {
    return false;
  }

  Atom::Rank id_count = 0;
  for (const auto&[element_name, count]:elements_counts)
  {
    double mass = elem_info::FindMass(element_name);
    for (int j = 0; j < count; ++j)
    {
      double position_X, position_Y, position_Z;
      if (!getline(ifs, line))
      {
        return false;
      }
      iss = std::istringstream(line);
      if (!(iss >> position_X >> position_Y >> position_Z))
      {
        return false;
      }
      atom_list_.emplace_back(id_count, mass, element_name,
                              position_X, position_Y, position_Z);
      element_list_set_[element_name].emplace_back(id_count);
      ++id_count;
    }
  }
  if (relOpt)
  {
    ConvertRelativeToAbsolute();
  } else
  {
    ConvertAbsoluteToRelative();
  }
  return true;
}

void Config::WriteConfig(const std::string &file_name) const
{
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "Number of particles = " << num_atoms_ << '\n';
  ofs << "A = " << scale_ << " Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << bravais_matrix_[kXDimension][kXDimension] << " A\n";
  ofs << "H0(1,2) = " << bravais_matrix_[kXDimension][kYDimension] << " A\n";
  ofs << "H0(1,3) = " << bravais_matrix_[kXDimension][kZDimension] << " A\n";
  ofs << "H0(2,1) = " << bravais_matrix_[kYDimension][kXDimension] << " A\n";
  ofs << "H0(2,2) = " << bravais_matrix_[kYDimension][kYDimension] << " A\n";
  ofs << "H0(2,3) = " << bravais_matrix_[kYDimension][kZDimension] << " A\n";
  ofs << "H0(3,1) = " << bravais_matrix_[kZDimension][kXDimension] << " A\n";
  ofs << "H0(3,2) = " << bravais_matrix_[kZDimension][kYDimension] << " A\n";
  ofs << "H0(3,3) = " << bravais_matrix_[kZDimension][kZDimension] << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom:atom_list_)
  {
    double mass = atom.GetMass();
    const std::string &type = atom.GetType();
    ofs << mass << '\n'
        << type << '\n'
        << atom.relative_position_[kXDimension] << " "
        << atom.relative_position_[kYDimension] << " "
        << atom.relative_position_[kZDimension] << '\n';
  }
}

void Config::WritePOSCAR(const std::string &file_name,
                         const bool &show_vacancy_option) const
{
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "#comment\n" << scale_ << '\n';
  ofs << bravais_matrix_[kXDimension][kXDimension] << " "
      << bravais_matrix_[kXDimension][kYDimension] << " "
      << bravais_matrix_[kXDimension][kZDimension] << '\n';
  ofs << bravais_matrix_[kYDimension][kXDimension] << " "
      << bravais_matrix_[kYDimension][kYDimension] << " "
      << bravais_matrix_[kYDimension][kZDimension] << '\n';
  ofs << bravais_matrix_[kZDimension][kXDimension] << " "
      << bravais_matrix_[kZDimension][kYDimension] << " "
      << bravais_matrix_[kZDimension][kZDimension] << '\n';
  std::ostringstream ele_oss, count_oss;
  for (const auto &[element, element_list]:element_list_set_)
  {
    if (!show_vacancy_option || element != "X")
    {
      ele_oss << element << " ";
      count_oss << element_list.size() << " ";
    }
  }
  ofs << ele_oss.str() << '\n' << count_oss.str() << '\n';
  ofs << "Direct\n";
  for (const auto &atom:atom_list_)
  {
    if (!show_vacancy_option || atom.GetType() != "X")
    {
      ofs << atom.relative_position_[kXDimension] << " "
          << atom.relative_position_[kYDimension] << " "
          << atom.relative_position_[kZDimension] << '\n';
    }
  }
}
void Config::GenerateUnitCell(const Matrix33 &bravais_matrix,
                              const std::vector<std::pair<std::string,
                                                          Vector3>> &type_position_list)
{
  Initialize();
  bravais_matrix_ = bravais_matrix;
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;
  for (const auto&[type, relative_position]:type_position_list)
  {
    atom_list_.emplace_back(num_atoms_,
                            elem_info::FindMass(type),
                            type,
                            relative_position);
    element_list_set_[type].emplace_back(num_atoms_++);
  }
  ConvertRelativeToAbsolute();
}
void Config::Duplicate(const std::array<int, kDimension> &factors)
{
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  bravais_matrix_ = {bravais_matrix_[kXDimension] * x_length,
                     bravais_matrix_[kYDimension] * y_length,
                     bravais_matrix_[kZDimension] * z_length};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  auto temp(std::move(atom_list_));
  num_atoms_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_set_.clear();
  neighbor_found_ = false;
  for (int k = 0; k < factors[kZDimension]; ++k)
  {
    for (int j = 0; j < factors[kYDimension]; ++j)
    {
      for (int i = 0; i < factors[kXDimension]; ++i)
      {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        for (const auto &atom:temp)
        {
          atom_list_.emplace_back(num_atoms_, atom.GetMass(), atom.GetType(),
                                  (x_reference + atom.relative_position_[kXDimension])
                                      / x_length,
                                  (y_reference + atom.relative_position_[kYDimension])
                                      / y_length,
                                  (z_reference + atom.relative_position_[kZDimension])
                                      / z_length);
          element_list_set_[atom.GetType()].emplace_back(num_atoms_++);
        }
      }
    }
  }

}
void Config::GenerateFCC(const double &lattice_constant_a,
                         const std::string &element,
                         const std::array<int, kDimension> &factors)
{
  Initialize();
  double mass = elem_info::FindMass(element);
  bravais_matrix_ = {{{lattice_constant_a * factors[kXDimension], 0, 0},
                      {0, lattice_constant_a * factors[kYDimension], 0},
                      {0, 0, lattice_constant_a * factors[kZDimension]}}};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k)
  {
    for (int j = 0; j < factors[kYDimension]; ++j)
    {
      for (int i = 0; i < factors[kXDimension]; ++i)
      {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 0.5) / x_length,
                                y_reference / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
void Config::GenerateBCC(const double &lattice_constant_a,
                         const std::string &element,
                         const std::array<int, kDimension> &factors)
{
  Initialize();
  double mass = elem_info::FindMass(element);
  bravais_matrix_ = {{{lattice_constant_a * factors[kXDimension], 0, 0},
                      {0, lattice_constant_a * factors[kYDimension], 0},
                      {0, 0, lattice_constant_a * factors[kZDimension]}}};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k)
  {
    for (int j = 0; j < factors[kYDimension]; ++j)
    {
      for (int i = 0; i < factors[kXDimension]; ++i)
      {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
void Config::GenerateHCP(const double &lattice_constant_a,
                         const double &lattice_constant_c,
                         const std::string &element,
                         const std::array<int, kDimension> &factors)
{
  Initialize();
  double mass = elem_info::FindMass(element);
  bravais_matrix_ = {{{lattice_constant_a * factors[kXDimension], 0, 0},
                      {-0.5 * lattice_constant_a * factors[kYDimension],
                       0.5 * sqrt(3) * lattice_constant_a * factors[kYDimension], 0},
                      {0, 0, lattice_constant_c * factors[kZDimension]}}};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;

  auto x_length = static_cast<double>(factors[kXDimension]);
  auto y_length = static_cast<double>(factors[kYDimension]);
  auto z_length = static_cast<double>(factors[kZDimension]);
  for (int k = 0; k < factors[kZDimension]; ++k)
  {
    for (int j = 0; j < factors[kYDimension]; ++j)
    {
      for (int i = 0; i < factors[kXDimension]; ++i)
      {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
        atom_list_.emplace_back(num_atoms_, mass, element,
                                (x_reference + 1.0 / 3.0) / x_length,
                                (y_reference + 2.0 / 3.0) / y_length,
                                (z_reference + 0.5) / z_length);
        element_list_set_[element].emplace_back(num_atoms_++);
      }
    }
  }
  ConvertRelativeToAbsolute();
}
const Matrix33 &Config::GetBravaisMatrix() const
{
  return bravais_matrix_;
}
const Matrix33 &Config::GetInverseBravaisMatrix() const
{
  return inverse_bravais_matrix_;
}
const Atom &Config::GetAtom(const Atom::Rank &index) const
{
  return atom_list_[index];
}
int Config::GetNumAtoms() const
{
  return num_atoms_;
}

}// namespace box
