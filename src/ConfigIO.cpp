#include"Config.h"

bool Config::ReadConfig(const std::string &file_name) {
  Initialize();
  std::ifstream ifs(file_name, std::ifstream::in);
  if (ifs.fail()) { return false; }
  std::string line;
  std::istringstream iss;
  double scale;
  Double3 first_bravais_vector{},
      second_bravais_vector{}, third_bravais_vector{};
  if (!getline(ifs, line)) { return false; }
  // "Number of particles = %i"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> num_atoms_)) { return false; }
  if (!getline(ifs, line)) { return false; }
  // A = 1.0 Angstrom (basic length-scale)
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> scale)) { return false; }
  box_.SetScale(scale);

  if (!getline(ifs, line)) { return false; }
  // "H0(1,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> first_bravais_vector[kXDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(1,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> first_bravais_vector[kYDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(1,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> first_bravais_vector[kZDim])) { return false; }
  box_.SetFirstBravaisVector(first_bravais_vector);

  if (!getline(ifs, line)) { return false; }
  // "H0(2,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> second_bravais_vector[kXDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(2,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> second_bravais_vector[kYDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(2,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> second_bravais_vector[kZDim])) { return false; }
  box_.SetSecondBravaisVector(second_bravais_vector);

  if (!getline(ifs, line)) { return false; }
  // "H0(3,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> third_bravais_vector[kXDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> third_bravais_vector[kYDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> third_bravais_vector[kZDim])) { return false; }
  box_.SetThirdBravaisVector(third_bravais_vector);

  if (!getline(ifs, line)) { return false; }
  // .NO_VELOCITY.
  if (!getline(ifs, line)) { return false; }
  // "entry_count = 3"
  for (int i = 0; i < num_atoms_; ++i) {
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
  double scale;
  Double3 first_bravais_vector{},
      second_bravais_vector{}, third_bravais_vector{};
  if (!getline(ifs, line)) { return false; }
  // #comment
  if (!getline(ifs, line)) { return false; }
  // scale factor, usually which is 1
  iss = std::istringstream(line);
  if (!(iss >> scale)) { return false; }
  box_.SetScale(scale);

  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> first_bravais_vector[kXDim] >> first_bravais_vector[kYDim]
            >> first_bravais_vector[kZDim])) { return false; }
  box_.SetFirstBravaisVector(first_bravais_vector);

  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> second_bravais_vector[kXDim] >> second_bravais_vector[kYDim]
            >> second_bravais_vector[kZDim])) { return false; }
  box_.SetSecondBravaisVector(second_bravais_vector);

  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> third_bravais_vector[kXDim] >> third_bravais_vector[kYDim]
            >> third_bravais_vector[kZDim])) { return false; }
  box_.SetThirdBravaisVector(third_bravais_vector);

  if (!getline(ifs, line)) { return false; }
  std::vector<std::string> elem_names;
  std::string elem;
  std::istringstream ele_iss(line);
  while (ele_iss >> elem) {
    elem_names.push_back(elem);
  }
  if (!getline(ifs, line)) { return false; }
  std::vector<int> elem_counts;
  int count;
  std::istringstream count_iss(line);
  while (count_iss >> count) { elem_counts.push_back(count); }
  num_atoms_ = accumulate(elem_counts.begin(), elem_counts.end(), 0);

  if (!getline(ifs, line)) { return false; }
  bool relOpt;
  if (line[0] == 'D' || line[0] == 'd') {
    relOpt = true;
  } else if (line[0] == 'C' || line[0] == 'c') {
    relOpt = false;
  } else {
    return false;
  }

  int id_count = 0;
  for (int i = 0; i < elem_names.size(); ++i) {
    double mass = elem_info::FindMass(elem_names[i]);
    for (int j = 0; j < elem_counts[i]; ++j) {
      double position_X, position_Y, position_Z;
      if (!getline(ifs, line)) { return false; }
      iss = std::istringstream(line);
      if (!(iss >> position_X >> position_Y >> position_Z)) { return false; }
      atom_list_.emplace_back(id_count++, mass, elem_names[i],
                              position_X, position_Y, position_Z);
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
  auto first_bravais_vector = box_.GetFirstBravaisVector();
  auto second_bravais_vector = box_.GetSecondBravaisVector();
  auto third_bravais_vector = box_.GetThirdBravaisVector();
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "Number of particles = " << num_atoms_ << "\n";
  ofs << "A = " << box_.GetScale() << " Angstrom (basic length-scale)\n";
  ofs << "H0(1,1) = " << first_bravais_vector[kXDim] << " A\n";
  ofs << "H0(1,2) = " << first_bravais_vector[kYDim] << " A\n";
  ofs << "H0(1,3) = " << first_bravais_vector[kZDim] << " A\n";
  ofs << "H0(2,1) = " << second_bravais_vector[kXDim] << " A\n";
  ofs << "H0(2,2) = " << second_bravais_vector[kYDim] << " A\n";
  ofs << "H0(2,3) = " << second_bravais_vector[kZDim] << " A\n";
  ofs << "H0(3,1) = " << third_bravais_vector[kXDim] << " A\n";
  ofs << "H0(3,2) = " << third_bravais_vector[kYDim] << " A\n";
  ofs << "H0(3,3) = " << third_bravais_vector[kZDim] << " A\n";
  ofs << ".NO_VELOCITY.\n";
  ofs << "entry_count = 3\n";
  for (const auto &atom : atom_list_) {
    double mass = atom.GetMass();
    const std::string &type = atom.GetType();
    auto relative_position = atom.GetRelativePosition();
    ofs << mass << "\n"
        << type << "\n"
        << relative_position[kXDim] << " "
        << relative_position[kYDim] << " "
        << relative_position[kZDim] << "\n";
  }
  ofs.close();
}
void Config::WritePOSCAR(const std::string &file_name,
                         const bool &show_vacancy_option) const {
  auto first_bravais_vector = box_.GetFirstBravaisVector();
  auto second_bravais_vector = box_.GetSecondBravaisVector();
  auto third_bravais_vector = box_.GetThirdBravaisVector();
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "#comment\n" << box_.GetScale() << "\n";
  ofs << first_bravais_vector[kXDim] << " "
      << first_bravais_vector[kYDim] << " "
      << first_bravais_vector[kZDim] << "\n";
  ofs << second_bravais_vector[kXDim] << " "
      << second_bravais_vector[kYDim] << " "
      << second_bravais_vector[kZDim] << "\n";
  ofs << third_bravais_vector[kXDim] << " "
      << third_bravais_vector[kYDim] << " "
      << third_bravais_vector[kZDim] << "\n";
  std::map<std::string, int> elem_counts;
  for (const auto &atm : atom_list_) {
    elem_counts[atm.GetType()]++;
  }
  std::ostringstream ele_oss, count_oss;
  for (const auto &elem_count:elem_counts) {
    if (!show_vacancy_option || elem_count.first != "Vac") {
      ele_oss << elem_count.first << " ";
      count_oss << elem_count.second << " ";
    }
  }
  ofs << ele_oss.str() << "\n" << count_oss.str() << "\n";
  ofs << "Direct\n";;
  for (const auto &atom : atom_list_) {
    if (!show_vacancy_option || atom.GetType() != "Vac") {
      auto relative_position = atom.GetRelativePosition();
      ofs << relative_position[kXDim] << " "
          << relative_position[kYDim] << " "
          << relative_position[kZDim] << "\n";
    }
  }
  ofs.close();
}

