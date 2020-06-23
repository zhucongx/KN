#include "ConfigIO.h"
#include <fstream>
#include <sstream>

namespace kn {

Config ConfigIO::ReadPOSCAR(const std::string &filename) {
  std::ifstream ifs(filename, std::ifstream::in);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // #comment
  double scale;
  ifs >> scale; // scale factor, usually which is 1.0
  Matrix33 basis;
  ifs >> basis;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line

  std::string buffer;
  getline(ifs, buffer);
  std::istringstream element_iss(buffer);
  getline(ifs, buffer);
  std::istringstream count_iss(buffer);

  std::string element;
  int count;
  int num_atoms = 0;
  std::vector<std::pair<std::string, int>> elements_counts;
  while (element_iss >> element && count_iss >> count) {
    elements_counts.emplace_back(element, count);
    num_atoms += count;
  }
  getline(ifs, buffer);
  bool relative_option;
  relative_option = buffer[0] == 'D' || buffer[0] == 'd';
  Config config(basis * scale, num_atoms);
  int id_count = 0;
  double position_X, position_Y, position_Z;
  for (const auto &[element_name, count] : elements_counts) {
    double mass = elem_info::FindMass(element_name);
    for (int j = 0; j < count; ++j) {
      ifs >> position_X >> position_Y >> position_Z;
      config.AppendAtom({id_count, mass, element_name,
                         position_X * scale, position_Y * scale, position_Z * scale});
      ++id_count;
    }
  }
  if (relative_option)
    config.ConvertRelativeToCartesian();
  else
    config.ConvertCartesianToRelative();

  return config;
}

Config ConfigIO::ReadConfig(const std::string &filename, bool update_neighbors) {
  std::ifstream ifs(filename, std::ifstream::in);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '='); // "Number of particles = %i"
  int num_atoms;
  ifs >> num_atoms;

  ifs.ignore(std::numeric_limits<std::streamsize>::max(),
             '='); // A = 1.0 Angstrom (basic length-scale)
  double scale;
  ifs >> scale;

  double basis_xx, basis_xy, basis_xz, basis_yx, basis_yy, basis_yz, basis_zx, basis_zy, basis_zz;
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
  Config config(Matrix33{{{basis_xx, basis_xy, basis_xz},
                          {basis_yx, basis_yy, basis_yz},
                          {basis_zx, basis_zy, basis_zz}}} * scale, num_atoms);

  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // .NO_VELOCITY.
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // "entry_count = 3"

  double mass, relative_position_X, relative_position_Y, relative_position_Z;
  int index;
  bool neighbor_found = false;
  for (int id = 0; id < num_atoms; ++id) {
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
      for (int i = 0; i < Al_const::kNumFirstNearestNeighbors; ++i) {
        ifs >> index;
        atom.first_nearest_neighbor_list_.push_back(index);
      }
      for (int i = 0; i < Al_const::kNumNearNeighbors; ++i) {
        ifs >> index;
        atom.near_neighbor_list_.push_back(index);
      }
      neighbor_found = true;
    }
    config.AppendAtom(atom);
  }
  config.ConvertRelativeToCartesian();
  config.neighbor_found_ = neighbor_found;
  if (update_neighbors)
    config.UpdateNeighbors(Al_const::kFirstNearestNeighborCutoff,
                           Al_const::kNearNeighborsCutoff);
  return config;
}

void ConfigIO::WritePOSCAR(const Config &config,
                           const std::string &filename,
                           bool show_vacancy_option) {
  std::ofstream ofs(filename, std::ofstream::out);
  ofs << "#comment\n1.0\n";
  ofs << config.GetBasis() << '\n';
  std::ostringstream ele_oss, count_oss;
  for (const auto &[element, element_list] : config.GetElementListMap()) {
    if (show_vacancy_option || element != "X") {
      ele_oss << element << ' ';
      count_oss << element_list.size() << ' ';
    }
  }
  ofs << ele_oss.str() << '\n' << count_oss.str() << '\n';
  ofs << "Direct\n";

  for (const auto &[element, element_list] : config.GetElementListMap()) {
    if (show_vacancy_option || element != "X") {
      for (auto index : element_list) {
        ofs << config.GetAtomList()[index].relative_position_ << '\n';
      }
    }
  }
}

void ConfigIO::WriteConfig(const Config &config, const std::string &filename, bool neighbors_info) {
  std::ofstream ofs(filename, std::ofstream::out);
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
    ofs << atom.mass_ << '\n'
        << atom.type_ << '\n'
        << atom.relative_position_;
    if (neighbors_info) {
      ofs << " #";
      for (auto neighbor_index : atom.first_nearest_neighbor_list_) {
        ofs << neighbor_index << ' ';
      }
      for (auto neighbor_index : atom.near_neighbor_list_) {
        ofs << neighbor_index << ' ';
      }
    }
    ofs << '\n';
  }
}
} // namespace kn
