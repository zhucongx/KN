//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _ATOM_H_
#define _ATOM_H_

#include <string>
#include <array>
#include <vector>

using std::vector;
using std::string;
using std::array;
// For FCC the first nearest neighbor is 12
#define firstNearestNbr 12
// 3D simulation
#define dimension 3

class Atom {
 private:
  // atom id which is an unique int for every atom indexed form 0
  int id;
  string type;
  // real position
  double position[dimension];
  // relative position in the box
  double relativePosition[dimension];
  //
  vector<int> nearNbrList;
  // First nearest neighbor list
  array<int, firstNearestNbr> firstNearestNbrList;
 public:
  Atom();
  Atom(int inId);
  Atom(int inId, double x, double y, double z);
  Atom(int inId, const string &inType, double x, double y, double z);
  ~Atom();
  bool operator<(const Atom &rhs) const;
  bool operator>(const Atom &rhs) const;
  bool operator<=(const Atom &rhs) const;
  bool operator>=(const Atom &rhs) const;
  int getId() const;
  void setId(int id);
  const string &getType() const;
  void setType(const string &type);
};

#endif //_ATOMCLASS_H_
