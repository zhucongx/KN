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

