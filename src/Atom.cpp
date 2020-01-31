//
// Created by Zhucong Xi on 1/31/20.
//

#include "Atom.h"

Atom::Atom() : id(0), type("X") {}
Atom::Atom(int inId) : id(inId), type("X") {}
Atom::Atom(int inId, double x, double y, double z) : id(inId), type("X") {
  relativePosition[0] = x;
  relativePosition[1] = y;
  relativePosition[2] = z;
  firstNearestNbrList.fill(-1);
}
Atom::Atom(int inId, const string &inType, double x, double y, double z)
    : id(inId), type(inType) {
  relativePosition[0] = x;
  relativePosition[1] = y;
  relativePosition[2] = z;
  firstNearestNbrList.fill(-1);
}
Atom::~Atom() {}
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
void Atom::setId(int id) {
  Atom::id = id;
}
const string &Atom::getType() const {
  return type;
}
void Atom::setType(const string &type) {
  Atom::type = type;
}

