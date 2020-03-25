#include"Config.h"

Config::Config() = default;
Config::~Config() = default;

void Config::clear(){
  num_atoms_ = 0;
  energy_ = 0;
  x_bravais_vector_.fill(0);
  y_bravais_vector_.fill(0);
  z_bravais_vector_.fill(0);
  atom_list_.clear();
}

bool Config::operator<(const Config &rhs) const {
  return energy_ < rhs.energy_;
}

void Config::ConvertRelativeToAbsolute() {
  for (auto &atom:atom_list_) {
    std::array<double, kDimension> absolute_position{};
    std::array<double, kDimension> relative_position =
        atom.GetRelativePosition();

    for (const auto &i : {kXDim, kYDim, kZDim}) {
      absolute_position[i] = relative_position[kXDim] * x_bravais_vector_[i]
          + relative_position[kYDim] * y_bravais_vector_[i]
          + relative_position[kZDim] * z_bravais_vector_[i];
    }
    atom.SetAbsolutePosition(absolute_position);
  }
}
void Config::ConvertAbsoluteToRelative() {
  arma::mat bm = {{x_bravais_vector_[kXDim], x_bravais_vector_[kYDim],
                   x_bravais_vector_[kZDim]},
                  {y_bravais_vector_[kXDim], y_bravais_vector_[kYDim],
                   y_bravais_vector_[kZDim]},
                  {z_bravais_vector_[kXDim], z_bravais_vector_[kYDim],
                   z_bravais_vector_[kZDim]}};
  for (auto &atom:atom_list_) {
    std::array<double, kDimension> absolute_position =
        atom.GetAbsolutePosition();
    std::array<double, kDimension> relative_position{};
    arma::vec b = {absolute_position[kXDim],
                   absolute_position[kYDim],
                   absolute_position[kZDim]};
    arma::vec x = solve(bm, b);
    for (const auto &i : {kXDim, kYDim, kZDim}) {
      relative_position[i] = x[i];
    }
    atom.SetRelativePosition(relative_position);
  }
}

bool Config::ReadConfig(const std::string &file_name) {
  clear();
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
  if (!getline(ifs, line)) { return false; }
  // "H0(1,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> x_bravais_vector_[kXDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(1,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> x_bravais_vector_[kYDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(1,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> x_bravais_vector_[kZDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(2,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> y_bravais_vector_[kXDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(2,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> y_bravais_vector_[kYDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // sscanf(line.c_str(), "H0(2,3) = %lf A", &y_bravais_vector_[kZDim]);
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> y_bravais_vector_[kZDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,1) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> z_bravais_vector_[kXDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,2) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> z_bravais_vector_[kYDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  // "H0(3,3) = %lf A"
  iss = std::istringstream(line);
  iss.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  if (!(iss >> z_bravais_vector_[kZDim])) { return false; }
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
  clear();
  std::ifstream ifs(file_name, std::ifstream::in);
  if (ifs.fail()) { return false; }
  std::string line;
  std::istringstream iss;
  if (!getline(ifs, line)) { return false; }
  // #comment
  if (!getline(ifs, line)) { return false; }
  // scale factor, usually which is 1
  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> x_bravais_vector_[kXDim] >> x_bravais_vector_[kYDim]
            >> x_bravais_vector_[kZDim])) { return false; }
  if (!getline(ifs, line)) { return false; }

  iss = std::istringstream(line);
  if (!(iss >> y_bravais_vector_[kXDim] >> y_bravais_vector_[kYDim]
            >> y_bravais_vector_[kZDim])) { return false; }
  if (!getline(ifs, line)) { return false; }
  iss = std::istringstream(line);
  if (!(iss >> z_bravais_vector_[kXDim] >> z_bravais_vector_[kYDim]
            >> z_bravais_vector_[kZDim])) { return false; }
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
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "Number of particles = " << num_atoms_ << std::endl;
  ofs << "A = 1.0 Angstrom (basic length-scale)" << std::endl;
  ofs << "H0(1,1) = " << x_bravais_vector_[0] << " A" << std::endl;
  ofs << "H0(1,2) = " << x_bravais_vector_[1] << " A" << std::endl;
  ofs << "H0(1,3) = " << x_bravais_vector_[2] << " A" << std::endl;
  ofs << "H0(2,1) = " << y_bravais_vector_[0] << " A" << std::endl;
  ofs << "H0(2,2) = " << y_bravais_vector_[1] << " A" << std::endl;
  ofs << "H0(2,3) = " << y_bravais_vector_[2] << " A" << std::endl;
  ofs << "H0(3,1) = " << z_bravais_vector_[0] << " A" << std::endl;
  ofs << "H0(3,2) = " << z_bravais_vector_[1] << " A" << std::endl;
  ofs << "H0(3,3) = " << z_bravais_vector_[2] << " A" << std::endl;
  ofs << ".NO_VELOCITY." << std::endl;
  ofs << "entry_count = 3" << std::endl;
  for (const auto &atom : atom_list_) {
    double mass = atom.GetMass();
    const std::string &type = atom.GetType();
    auto relative_position = atom.GetRelativePosition();
    ofs << ((mass > 0) ? mass : elem_info::FindMass(type)) << std::endl
        << type << std::endl
        << relative_position[kXDim] << " " << relative_position[kYDim] << " "
        << relative_position[kZDim] << std::endl;
  }
  ofs.close();
}
void Config::WritePOSCAR(const std::string &file_name,
                         const bool &show_vacancy_option) const {
  std::ofstream ofs(file_name, std::ofstream::out);
  ofs << "#comment" << std::endl << "1.00000" << std::endl;
  ofs << x_bravais_vector_[kXDim] << " "
      << x_bravais_vector_[kYDim] << " "
      << x_bravais_vector_[kZDim] << std::endl;
  ofs << y_bravais_vector_[kXDim] << " "
      << y_bravais_vector_[kYDim] << " "
      << y_bravais_vector_[kZDim] << std::endl;
  ofs << z_bravais_vector_[kXDim] << " "
      << z_bravais_vector_[kYDim] << " "
      << z_bravais_vector_[kZDim] << std::endl;
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
  ofs << ele_oss.str() << std::endl << count_oss.str() << std::endl;
  ofs << "Direct" << std::endl;
  for (const auto &atom : atom_list_) {
    if (!show_vacancy_option || atom.GetType() != "Vac") {
      auto relative_position = atom.GetRelativePosition();
      ofs << relative_position[kXDim] << " " << relative_position[kYDim] << " "
          << relative_position[kZDim] << std::endl;
    }
  }
  ofs.close();
}

