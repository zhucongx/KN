//
//CreatedbyZhucongXion1/31/20.
//

#include"Configuration.h"

Configuration::Configuration() : numAtoms(0), numTypes(0), energy(0.0) {}
Configuration::~Configuration() {}
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
void Configuration::readLammpsData(const string &fileName) {
}
bool Configuration::readConfig(const string &fileName) {
  ifstream ifs(fileName, ifstream::in);
  if (ifs.fail()) {
    std::cout << "Error when open file." << std::endl;
    return false;
  }
  string buff;
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "Number of particles = %i", &numAtoms);
  std::getline(ifs, buff); // A = 1.0 Angstrom (basic length-scale)
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,1) = %lf A", &bvx[X]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,2) = %lf A", &bvx[Y]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(1,3) = %lf A", &bvx[Z]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,1) = %lf A", &bvy[X]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,2) = %lf A", &bvy[Y]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(2,3) = %lf A", &bvy[Z]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,1) = %lf A", &bvz[X]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,2) = %lf A", &bvz[Y]);
  std::getline(ifs, buff);
  sscanf(buff.c_str(), "H0(3,3) = %lf A", &bvz[Z]);
  std::getline(ifs, buff); // .NO_VELOCITY.
  std::getline(ifs, buff);
  int entry = 3;
  sscanf(buff.c_str(), "entry_count = %i", &entry);
}

