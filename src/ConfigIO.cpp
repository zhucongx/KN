#include"Config.h"

namespace box {

bool Config::ReadConfig(const std::string &file_name) {
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);
  if (ifs.fail()) { return false; }
  std::string line;
  std::istringstream iss;
  if (!getline(ifs, line)) { return false; }
  // "Number of particles = %i"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> num_atoms_)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // A = 1.0 Angstrom (basic length-scale)
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> scale_)) { return false; }

  if (!getline(ifs, line)) { return false; }
  // "H0(1,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> first_bravais_vector_.x)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(1,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> first_bravais_vector_.y)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(1,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> first_bravais_vector_.z)) { return false; }

  if (!getline(ifs, line)) { return false; }
  // "H0(2,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> second_bravais_vector_.x)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(2,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> second_bravais_vector_.y)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(2,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> second_bravais_vector_.z)) { return false; }

  if (!getline(ifs, line)) { return false; }
  // "H0(3,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> third_bravais_vector_.x)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> third_bravais_vector_.y)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> third_bravais_vector_.z)) { return false; }

  if (!getline(ifs, line)) { return false; }
  // .NO_VELOCITY.
  if (!getline(ifs, line)) { return false; }
  // "entry_count = 3"
  for (Atom::Rank i = 0; i < num_atoms_; ++i) {
    double mass, relative_position_X, relative_position_Y, relative_position_Z;
    std::string type;
    if (!getline(ifs, line)) { return false; }
    iss = std::istringstream(line);
    if (!(iss >> mass)) { return false; }
    if (!getline(ifs, line)) { return false; }
    type = line;
    if (!getline(ifs, line)) { return false; }
    iss = std::istringstream(line);
    if (!(iss >> relative_position_X >> relative_position_Y
              >> relative_position_Z)) { return false; }
    atom_list_.emplace_back(i, mass, type,
                            relative_position_X,
                            relative_position_Y,
                            relative_position_Z);
    element_list_set_[type].emplace_back(i);
  }
  ifs.close();
  ConvertRelativeToAbsolute();
  return true;
}
bool Config::ReadPOSCAR(const std::string &file_name) {
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);
  if (ifs.fail()) { return false; }
  std::string line;
  std::istringstream iss;
  if (!getline(ifs, line)) { return false; }
  // #comment
  if (!getline(ifs, line)) { return false; }
  // scale factor, usually which is 1
  iss = std::istringstream(line);
  if (!(iss >> scale_)) { return false; }

  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> first_bravais_vector_.x >> first_bravais_vector_.y
            >> first_bravais_vector_.z)) { return false; }

  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> second_bravais_vector_.x >> second_bravais_vector_.y
            >> second_bravais_vector_.z)) { return false; }

  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> third_bravais_vector_.x >> third_bravais_vector_.y
            >> third_bravais_vector_.z)) { return false; }

  if (!getline(ifs, line)) { return false; }
  std::string element;
  std::istringstream name_iss(line);
  if (!getline(ifs, line)) { return false; }
  int count;
  std::istringstream count_iss(line);

  std::vector<std::pair<std::string, int>> elements_counts;
  while (name_iss >> element && count_iss >> count) {
    num_atoms_ += count;
    elements_counts.emplace_back(element, count);
  }

  if (!getline(ifs, line)) { return false; }
  bool relOpt;
  if (line[0] == 'D' || line[0] == 'd') {
    relOpt = true;
  } else if (line[0] == 'C' || line[0] == 'c') {
    relOpt = false;
  } else {
    return false;
  }

  Atom::Rank id_count = 0;
  for (const auto&[element_name, count]:elements_counts) {
    double mass = elem_info::FindMass(element_name);
    for (int j = 0; j < count; ++j) {
      double position_X, position_Y, position_Z;
      if (!getline(ifs, line)) { return false; }
      iss = std::istringstream(line);
      if (!(iss >> position_X >> position_Y >> position_Z)) { return false; }
      atom_list_.emplace_back(id_count, mass, element_name,
                              position_X, position_Y, position_Z);
      element_list_set_[element_name].emplace_back(id_count);
      ++id_count;
    }
  }
  if (relOpt) {
    ConvertRelativeToAbsolute();
  } else {
    ConvertAbsoluteToRelative();
  }
  ifs.close();
  return true;
}
void Config::WriteConfig(const std::string &file_name) const {
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "Number of particles = " << num_atoms_ << "\n";
  ofs << "A = " << scale_ << " Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << first_bravais_vector_.x << " A\n";
  ofs << "H0(1,2) = " << first_bravais_vector_.y << " A\n";
  ofs << "H0(1,3) = " << first_bravais_vector_.z << " A\n";
  ofs << "H0(2,1) = " << second_bravais_vector_.x << " A\n";
  ofs << "H0(2,2) = " << second_bravais_vector_.y << " A\n";
  ofs << "H0(2,3) = " << second_bravais_vector_.z << " A\n";
  ofs << "H0(3,1) = " << third_bravais_vector_.x << " A\n";
  ofs << "H0(3,2) = " << third_bravais_vector_.y << " A\n";
  ofs << "H0(3,3) = " << third_bravais_vector_.z << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom:atom_list_) {
    double mass = atom.GetMass();
    const std::string &type = atom.GetType();
    ofs << mass << "\n"
        << type << "\n"
        << atom.relative_position_.x << " "
        << atom.relative_position_.y << " "
        << atom.relative_position_.z << "\n";
  }
  ofs.close();
}
void Config::WritePOSCAR(const std::string &file_name,
                         const bool &show_vacancy_option) const {
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "#comment\n" << scale_ << "\n";
  ofs << first_bravais_vector_.x << " "
      << first_bravais_vector_.y << " "
      << first_bravais_vector_.z << "\n";
  ofs << second_bravais_vector_.x << " "
      << second_bravais_vector_.y << " "
      << second_bravais_vector_.z << "\n";
  ofs << third_bravais_vector_.x << " "
      << third_bravais_vector_.y << " "
      << third_bravais_vector_.z << "\n";
  std::ostringstream ele_oss, count_oss;
  for (const auto &[element, element_list]:element_list_set_) {
    if (!show_vacancy_option || element != "Vac") {
      ele_oss << element << " ";
      count_oss << element_list.size() << " ";
    }
  }
  ofs << ele_oss.str() << "\n" << count_oss.str() << "\n";
  ofs << "Direct\n";
  for (const auto &atom:atom_list_) {
    if (!show_vacancy_option || atom.GetType() != "Vac") {
      ofs << atom.relative_position_.x << " "
          << atom.relative_position_.y << " "
          << atom.relative_position_.z << "\n";
    }
  }
  ofs.close();
}

}// namespace box
