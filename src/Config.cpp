#include"Config.h"

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
  bravais_matrix_ = {{0, 0, 0},
                     {0, 0, 0},
                     {0, 0, 0}};
  num_atoms_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_set_.clear();
  neighbor_found_ = false;
}

bool Config::IsCubic() const
{
  return bravais_matrix_.row1.x == bravais_matrix_.row2.y &&
      bravais_matrix_.row2.y == bravais_matrix_.row3.z &&
      bravais_matrix_.row1.y == 0 &&
      bravais_matrix_.row1.z == 0 &&
      bravais_matrix_.row2.x == 0 &&
      bravais_matrix_.row2.z == 0 &&
      bravais_matrix_.row3.x == 0 &&
      bravais_matrix_.row3.y == 0;
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
    add_displacement(atom.cartesian_position_.x);
    add_displacement(atom.cartesian_position_.y);
    add_displacement(atom.cartesian_position_.z);
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
      Vector3<double> absolute_distance_vector =
          GetRelativeDistanceVector(*it1, *it2) * bravais_matrix_;
      if (absolute_distance_vector.x > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector.y > second_r_cutoff_square)
        continue;
      if (absolute_distance_vector.z > second_r_cutoff_square)
        continue;
      double absolute_distance_square =
          InnerProduct(absolute_distance_vector);
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
    atom.relative_position_.x -= floor(atom.relative_position_.x);
    atom.relative_position_.y -= floor(atom.relative_position_.y);
    atom.relative_position_.z -= floor(atom.relative_position_.z);

    atom.cartesian_position_ = atom.relative_position_ * bravais_matrix_;
  }
}

// void Config::WrapAbsolutePosition() {
//   ConvertAbsoluteToRelative();
//   WrapRelativePosition();
// }

// void Config::ShiftAtomToCentral(const Atom::Rank &id) {
//   Vector3<double> criterionDistance =
//       atom_list_[id].relative_position_ - Vector3<double>{0.5, 0.5, 0.5};
//   for (auto &atom:atom_list_) {
//     atom.relative_position_ = atom.relative_position_ - criterionDistance;
//   }
//   WrapRelativePosition();
// }
// for better performance, shouldn't call Wrap function
void Config::MoveRelativeDistance(const Vector3<double> &distance_vector)
{
  for (auto &atom:atom_list_)
  {
    atom.relative_position_ += distance_vector;

    atom.relative_position_ -= Floor(atom.relative_position_);

    atom.cartesian_position_ = atom.relative_position_ * bravais_matrix_;
  }
}
void Config::MoveOneAtomRelativeDistance(const Atom::Rank &index,
                                         const Vector3<double> &distance_vector)
{
  atom_list_[index].relative_position_ += distance_vector;
  atom_list_[index].relative_position_ -=
      Floor(atom_list_[index].relative_position_);

  atom_list_[index].cartesian_position_ =
      atom_list_[index].relative_position_ * bravais_matrix_;
}
// void Config::MoveAbsoluteDistance(const Vector3<double> &distance_vector) {
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
  if (!getline(ifs, line))
  {
    return false;
  }
  // "Number of particles = %i"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> num_atoms_))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // A = 1.0 Angstrom (basic length-scale)
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> scale_))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(1,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row1.x))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(1,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row1.y))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(1,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row1.z))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(2,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row2.x))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(2,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row2.y))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(2,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row2.z))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(3,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row3.x))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(3,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row3.y))
  {
    return false;
  }
  if (!getline(ifs, line))
  {
    return false;
  }
  // "H0(3,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> bravais_matrix_.row3.z))
  {
    return false;
  }
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);

  if (!getline(ifs, line))
  {
    return false;
  }
  // .NO_VELOCITY.
  if (!getline(ifs, line))
  {
    return false;
  }
  // "entry_count = 3"
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
  if (!getline(ifs, line))
  {
    return false;
  }
  // #comment
  if (!getline(ifs, line))
  {
    return false;
  }
  // scale factor, usually which is 1
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
  if (!(iss >> bravais_matrix_.row1.x >> bravais_matrix_.row1.y
            >> bravais_matrix_.row1.z))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  iss = std::istringstream(line);
  if (!(iss >> bravais_matrix_.row2.x >> bravais_matrix_.row2.y
            >> bravais_matrix_.row2.z))
  {
    return false;
  }

  if (!getline(ifs, line))
  {
    return false;
  }
  iss = std::istringstream(line);
  if (!(iss >> bravais_matrix_.row3.x >> bravais_matrix_.row3.y
            >> bravais_matrix_.row3.z))
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
  ofs << "H0(1,1) = " << bravais_matrix_.row1.x << " A\n";
  ofs << "H0(1,2) = " << bravais_matrix_.row1.y << " A\n";
  ofs << "H0(1,3) = " << bravais_matrix_.row1.z << " A\n";
  ofs << "H0(2,1) = " << bravais_matrix_.row2.x << " A\n";
  ofs << "H0(2,2) = " << bravais_matrix_.row2.y << " A\n";
  ofs << "H0(2,3) = " << bravais_matrix_.row2.z << " A\n";
  ofs << "H0(3,1) = " << bravais_matrix_.row3.x << " A\n";
  ofs << "H0(3,2) = " << bravais_matrix_.row3.y << " A\n";
  ofs << "H0(3,3) = " << bravais_matrix_.row3.z << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom:atom_list_)
  {
    double mass = atom.GetMass();
    const std::string &type = atom.GetType();
    ofs << mass << '\n'
        << type << '\n'
        << atom.relative_position_.x << " "
        << atom.relative_position_.y << " "
        << atom.relative_position_.z << '\n';
  }
}

void Config::WritePOSCAR(const std::string &file_name,
                         const bool &show_vacancy_option) const
{
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "#comment\n" << scale_ << '\n';
  ofs << bravais_matrix_.row1.x << " "
      << bravais_matrix_.row1.y << " "
      << bravais_matrix_.row1.z << '\n';
  ofs << bravais_matrix_.row2.x << " "
      << bravais_matrix_.row2.y << " "
      << bravais_matrix_.row2.z << '\n';
  ofs << bravais_matrix_.row3.x << " "
      << bravais_matrix_.row3.y << " "
      << bravais_matrix_.row3.z << '\n';
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
      ofs << atom.relative_position_.x << " "
          << atom.relative_position_.y << " "
          << atom.relative_position_.z << '\n';
    }
  }
}
void Config::GenerateUnitCell(const Matrix33<double> &bravais_matrix,
                              const std::vector<std::pair<std::string,
                                                          Vector3<double>>> &type_position_list)
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
void Config::Duplicate(const Vector3<int> &factors)
{
  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  bravais_matrix_ = {bravais_matrix_.row1 * x_length,
                     bravais_matrix_.row2 * y_length,
                     bravais_matrix_.row3 * z_length};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  auto temp(std::move(atom_list_));
  num_atoms_ = 0;
  energy_ = 0.0;
  scale_ = 1.0;
  atom_list_.clear();
  element_list_set_.clear();
  neighbor_found_ = false;
  for (int k = 0; k < factors.z; ++k)
  {
    for (int j = 0; j < factors.y; ++j)
    {
      for (int i = 0; i < factors.x; ++i)
      {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        for (const auto &atom:temp)
        {
          atom_list_.emplace_back(num_atoms_, atom.GetMass(), atom.GetType(),
                                  (x_reference + atom.relative_position_.x)
                                      / x_length,
                                  (y_reference + atom.relative_position_.y)
                                      / y_length,
                                  (z_reference + atom.relative_position_.z)
                                      / z_length);
          element_list_set_[atom.GetType()].emplace_back(num_atoms_++);
        }
      }
    }
  }

}
void Config::GenerateFCC(const double &lattice_constant_a,
                         const std::string &element,
                         const Vector3<int> &factors)
{
  Initialize();
  double mass = elem_info::FindMass(element);
  bravais_matrix_ = {{lattice_constant_a * factors.x, 0, 0},
                     {0, lattice_constant_a * factors.y, 0},
                     {0, 0, lattice_constant_a * factors.z}};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  for (int k = 0; k < factors.z; ++k)
  {
    for (int j = 0; j < factors.y; ++j)
    {
      for (int i = 0; i < factors.x; ++i)
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
                         const Vector3<int> &factors)
{
  Initialize();
  double mass = elem_info::FindMass(element);
  bravais_matrix_ = {{lattice_constant_a * factors.x, 0, 0},
                     {0, lattice_constant_a * factors.y, 0},
                     {0, 0, lattice_constant_a * factors.z}};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  for (int k = 0; k < factors.z; ++k)
  {
    for (int j = 0; j < factors.y; ++j)
    {
      for (int i = 0; i < factors.x; ++i)
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
                         const Vector3<int> &factors)
{
  Initialize();
  double mass = elem_info::FindMass(element);
  bravais_matrix_ = {{lattice_constant_a * factors.x, 0, 0},
                     {-0.5 * lattice_constant_a * factors.y,
                      0.5 * sqrt(3) * lattice_constant_a * factors.y, 0},
                     {0, 0, lattice_constant_c * factors.z}};
  inverse_bravais_matrix_ = InverseMatrix33(bravais_matrix_);
  scale_ = 1.0;
  num_atoms_ = 0;

  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  for (int k = 0; k < factors.z; ++k)
  {
    for (int j = 0; j < factors.y; ++j)
    {
      for (int i = 0; i < factors.x; ++i)
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
const Matrix33<double> &Config::GetBravaisMatrix() const
{
  return bravais_matrix_;
}
const Matrix33<double> &Config::GetInverseBravaisMatrix() const
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
