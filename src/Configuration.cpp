//
//CreatedbyZhucongXion1/31/20.
//

#include"Configuration.h"

Configuration::Configuration() : numAtoms(0), numTypes(0), energy(0.0) {}
Configuration::~Configuration() = default;
bool Configuration::operator<(const Configuration &rhs) const {
  return energy < rhs.energy;
}
bool Configuration::operator>(const Configuration &rhs) const {
  return rhs < *this;
}
bool Configuration::operator<=(const Configuration &rhs) const {
  return !(rhs < *this);
}
bool Configuration::operator>=(const Configuration &rhs) const {
  return !(*this < rhs);
}
void Configuration::cnvPrl2Pst() {
  for (auto &atm:atoms) {
    atm.cnvPrl2Pst(bvx, bvy, bvz);
  }
}
void Configuration::cnvPst2Prl() {
  arma::mat bm = {{bvx[0], bvx[1], bvx[2]},
                  {bvy[0], bvy[1], bvy[2]},
                  {bvz[0], bvz[1], bvz[2]}};

  for (auto &atm:atoms) {
    atm.cnvPst2Prl(bm);
  }
}
bool Configuration::readLammpsData(const std::string &fileName) {}
bool Configuration::readConfig(const std::string &fileName) {
  std::ifstream ifs(fileName, std::ifstream::in);
  if (ifs.fail()) { return false; }
  std::string line;
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "Number of particles = %i", &numAtoms);
  if (!getline(ifs, line)) { return false; }
  // A = 1.0 Angstrom (basic length-scale)
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(1,1) = %lf A", &bvx[X]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(1,2) = %lf A", &bvx[Y]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(1,3) = %lf A", &bvx[Z]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(2,1) = %lf A", &bvy[X]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(2,2) = %lf A", &bvy[Y]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(2,3) = %lf A", &bvy[Z]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(3,1) = %lf A", &bvz[X]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(3,2) = %lf A", &bvz[Y]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "H0(3,3) = %lf A", &bvz[Z]);
  if (!getline(ifs, line)) { return false; }
  // .NO_VELOCITY.
  if (!getline(ifs, line)) { return false; }
  int entry = 3;
  sscanf(line.c_str(), "entry_count = %i", &entry);

  for (int i = 0; i < numAtoms; i++) {
    Atom atm;
    if (!atm.readConfig(ifs)) { return false; }
    atoms.push_back(atm);
  }
  ifs.close();
  cnvPrl2Pst();
  return true;
}
bool Configuration::readPOSCAR(const std::string &fileName) {
  std::ifstream ifs(fileName, std::ifstream::in);
  if (ifs.fail()) { return false; }
  std::string line;
  if (!getline(ifs, line)) { return false; }
  // #comment
  if (!getline(ifs, line)) { return false; }
  // scale factor, usually which is 1
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf %lf %lf", &bvx[X], &bvx[Y], &bvx[Z]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf %lf %lf", &bvy[X], &bvy[Y], &bvy[Z]);
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf %lf %lf", &bvz[X], &bvz[Y], &bvz[Z]);

  if (!getline(ifs, line)) { return false; }
  std::vector<std::string> elemNames;
  std::string elem;
  std::istringstream eleIss(line);
  while (eleIss >> elem) {
    elemNames.push_back(elem);
  }

  if (!getline(ifs, line)) { return false; }
  std::vector<int> elemCounts;
  int count;
  std::istringstream countIss(line);
  while (countIss >> count) { elemCounts.push_back(count); }
  numAtoms = accumulate(elemCounts.begin(), elemCounts.end(), 0);

  if (!getline(ifs, line)) { return false; }
  bool relOpt;
  if (line[0] == 'D' || line[0] == 'd') { relOpt = true; }
  else if (line[0] == 'C' || line[0] == 'c') { relOpt = false; }
  else { return false; }

  int idCount = 0;
  for (int i = 0; i < elemCounts.size(); i++) {
    for (int j = 0; j < elemCounts[i]; j++) {
      Atom atm(idCount, elemNames[i]);
      idCount++;
      atm.readPOSCAR(ifs, relOpt);
      atoms.push_back(atm);
    }
  }
  if (relOpt) { cnvPrl2Pst(); } else { cnvPst2Prl(); }
  ifs.close();
  return true;
}
void Configuration::writeConfig(const std::string &fileName) const {
  std::ofstream ofs(fileName, std::ofstream::out);
  ofs << "Number of particles = " << numAtoms << std::endl;
  ofs << "A = 1.0 Angstrom (basic length-scale)" << std::endl;
  ofs << "H0(1,1) = " << bvx[X] << " A" << std::endl;
  ofs << "H0(1,2) = " << bvx[Y] << " A" << std::endl;
  ofs << "H0(1,3) = " << bvx[Z] << " A" << std::endl;
  ofs << "H0(2,1) = " << bvy[X] << " A" << std::endl;
  ofs << "H0(2,2) = " << bvy[Y] << " A" << std::endl;
  ofs << "H0(2,3) = " << bvy[Z] << " A" << std::endl;
  ofs << "H0(3,1) = " << bvz[X] << " A" << std::endl;
  ofs << "H0(3,2) = " << bvz[Y] << " A" << std::endl;
  ofs << "H0(3,3) = " << bvz[Z] << " A" << std::endl;
  ofs << ".NO_VELOCITY." << std::endl;
  ofs << "entry_count = 3" << std::endl;
  for (const auto &atm : atoms) {
    atm.writeConfig(ofs);
  }
  ofs.close();
}
void Configuration::writePOSCAR(const std::string &fileName,
                                const bool &vacOption) const {
  std::ofstream ofs(fileName, std::ofstream::out);
  ofs << "#comment" << std::endl << "1.00000" << std::endl;
  ofs << bvx[X] << " " << bvx[Y] << " " << bvx[Z] << std::endl;
  ofs << bvy[X] << " " << bvy[Y] << " " << bvy[Z] << std::endl;
  ofs << bvz[X] << " " << bvz[Y] << " " << bvz[Z] << std::endl;

  std::map<std::string, int> elemCounts;
  for (const auto &atm : atoms) {
    elemCounts[atm.getType()]++;
  }

  std::ostringstream eleOss, countOss;
  for (const auto &elemCount:elemCounts) {
    if (!vacOption || elemCount.first != "X") {
      eleOss << elemCount.first << " ";
      countOss << elemCount.second << " ";
    }
  }
  ofs << eleOss.str() << std::endl << countOss.str() << std::endl;
  ofs << "Direct"<<std::endl;
  for (const auto &atm : atoms) {
    if (!vacOption || atm.getType() != "X") {
      atm.writePrl(ofs);
    }
  }
  ofs.close();
}





