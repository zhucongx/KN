//
// Created by Zhucong Xi on 1/31/20.
//

#include "Atom.h"

Atom::Atom() : id(0), type("X") {}
Atom::Atom(int i) : id(i), type("X") {}
Atom::Atom(int i, std::string tp) : id(i), type(std::move(tp)) {}
// Atom::Atom(int i, double x, double y, double z) : id(i), type("X") {
//   prl[0] = x;
//   prl[1] = y;
//   prl[2] = z;
//   firstNearestNbrList.fill(-1);
// }
Atom::Atom(int i, double m, std::string tp, double x, double y, double z)
    : id(i), mass(m), type(std::move(tp)) {
  prl[0] = x;
  prl[1] = y;
  prl[2] = z;
  pst[0] = x;
  pst[1] = y;
  pst[2] = z;
  firstNearestNbrList.fill(-1);
}
Atom::~Atom() = default;
bool Atom::operator<(const Atom &rhs) const {
  return type < rhs.type;
}
bool Atom::operator>(const Atom &rhs) const {
  return rhs < *this;
}
bool Atom::operator<=(const Atom &rhs) const {
  return !(rhs < *this);
}
bool Atom::operator>=(const Atom &rhs) const {
  return !(*this < rhs);
}

int Atom::getId() const {
  return id;
}
void Atom::setId(int i) {
  id = i;
}
const std::string &Atom::getType() const {
  return type;
}
void Atom::setType(const std::string &tp) {
  type = tp;
}
double Atom::findMass() const {
  /*this "it" is indeed a iterator*/
  auto it = std::find(elementList.begin(), elementList.end(), type);
  int index = std::distance(elementList.begin(), it);
  return massList[index];
}

void Atom::cnvPrl2Pst(const std::array<double, 3> &bvx,
                      const std::array<double, 3> &bvy,
                      const std::array<double, 3> &bvz) {
  for (const int i : {X, Y, Z}) {
    pst[i] = prl[X] * bvx[i] + prl[Y] * bvy[i] + prl[Z] * bvz[i];
  }
}

void Atom::cnvPst2Prl(const arma::mat &bm) {
  arma::vec b = {pst[X],
                 pst[Y],
                 pst[Z]};
  arma::vec x = solve(bm, b);
  prl[X] = x[X];
  prl[Y] = x[Y];
  prl[Z] = x[Z];
}

bool Atom::readConfig(std::ifstream &ifs) {
  std::string line;
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf", &mass);
  if (!getline(ifs, line)) { return false; }
  type = line;
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf %lf %lf", &prl[X],
         &prl[Y], &prl[Z]);
  return true;
}
void Atom::writeConfig(std::ofstream &ofs) const {
  ofs << ((mass > 0) ? mass : findMass()) << std::endl << type << std::endl
      << prl[X] << " " << prl[Y] << " " << prl[Z] << std::endl;
}
bool Atom::readPOSCARDirect(std::ifstream &ifs) {
  std::string line;
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf %lf %lf", &prl[X], &prl[Y], &prl[Z]);
  return true;
}
bool Atom::readPOSCARCartesian(std::ifstream &ifs) {
  std::string line;
  if (!getline(ifs, line)) { return false; }
  sscanf(line.c_str(), "%lf %lf %lf", &pst[X], &pst[Y], &pst[Z]);
  return true;
}

