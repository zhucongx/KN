//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _ATOM_H_
#define _ATOM_H_

#include <string>
#include <fstream>
#include <array>
#include <vector>
#include "armadillo"

// For FCC the first nearest neighbor is 12
#define firstNearestNbr 12
// 3D simulation
#define dimension 3
enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };

class Atom {
 private:
  // atom id which is an unique int for every atom indexed form 0
  int id;
  double mass{};
  std::string type;
  // real position
  double pst[dimension]{};
  // relative position in the box
  double prl[dimension]{};
  //
  std::vector<int> nearNbrList;
  // First nearest neighbor list
  std::array<int, firstNearestNbr> firstNearestNbrList{};
 public:
  Atom();
  Atom(int i);
  Atom(int i, std::string tp);
  // Atom(int i, double x, double y, double z);
  Atom(int i, double m, std::string tp, double x, double y, double z);
  virtual ~Atom();
  bool operator<(const Atom &rhs) const;
  bool operator>(const Atom &rhs) const;
  bool operator<=(const Atom &rhs) const;
  bool operator>=(const Atom &rhs) const;

  int getId() const;
  void setId(int i);
  const std::string &getType() const;
  void setType(const std::string &tp);
  double findMass() const;

  void cnvPrl2Pst(const std::array<double, 3> &bvx,
                  const std::array<double, 3> &bvy,
                  const std::array<double, 3> &bvz);
  void cnvPst2Prl(const arma::mat&);

  bool readConfig(std::ifstream &ifs);
  void writeConfig(std::ofstream &ofs) const;
  bool readPOSCAR(std::ifstream &ifs, const bool&realOption);
  void writePrl(std::ofstream &ofs) const;
};

#include "Elem.inl"

#endif //_ATOMCLASS_H_
