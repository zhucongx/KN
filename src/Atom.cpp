//
// Created by Zhucong Xi on 1/31/20.
//
#include "Atom.h"

Atom::Atom() : id(0), type("X") {}
Atom::Atom(int i) : id(i), type("X") {}
Atom::Atom(int i, double m, std::string tp)
    : id(i), mass(m), type(std::move(tp)) {}
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
  mass = getMass();
}
double Atom::getMass() const {
  return ElemInfo::findMass(type);
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
  sscanf(line.c_str(), "%lf %lf %lf", &prl[X], &prl[Y], &prl[Z]);
  return true;
}
// Read from POSCAR file to atom. If Direct meaning relative position,
// relativeOption should be true. If Cartesian meaning real position,
// relativeOption should be false.
bool Atom::readPOSCAR(std::ifstream &ifs, const bool &relativeOption) {
  std::string line;
  if (!getline(ifs, line)) { return false; }
  if (relativeOption) {
    sscanf(line.c_str(), "%lf %lf %lf", &prl[X], &prl[Y], &prl[Z]);
  } else {
    sscanf(line.c_str(), "%lf %lf %lf", &pst[X], &pst[Y], &pst[Z]);
  }
  return true;
}
void Atom::writeConfig(std::ofstream &ofs) const {
  ofs << ((mass > 0) ? mass : getMass()) << std::endl << type << std::endl;
  writePrl(ofs);
}
void Atom::writePrl(std::ofstream &ofs) const {
  ofs << prl[X] << " " << prl[Y] << " " << prl[Z] << std::endl;
}