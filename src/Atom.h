//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef KN_INCLUDE_ATOM_H_
#define KN_INCLUDE_ATOM_H_

#include <fstream>
#include <array>

#include "armadillo"

#include "ElemInfo.h"

// typedef int Rank;
// For FCC the first nearest neighbor is 12
const int kNumFirstNearestNeighbors = 12;
// For FCC the first and second nearest neighbor is 18
const int kNumNearNeighbors = 18;
// 3D simulation
const int kDimension = 3;
enum { kXDim = 0, kYDim = 1, kZDim = 2 };

class Atom {
 public:
  Atom(int id);
  Atom(int id, double mass, std::string type);
  // Atom(int i, double x, double y, double z);
  Atom(int i, double m, std::string tp, double x, double y, double z);
  virtual ~Atom();
  void SetId(int id);
  [[nodiscard]] int GetId() const;
  void SetMass(double mass);
  [[nodiscard]] double GetMass() const;
  void SetType(const std::string &type);
  [[nodiscard]] const std::string &GetType() const;
  void SetAbsolutePosition(const std::array<double,
                                            kDimension> &absolute_position);
  [[nodiscard]] const std::array<double,
                                 kDimension> &GetAbsolutePosition() const;
  void SetRelativePosition(const std::array<double,
                                            kDimension> &relative_position);
  [[nodiscard]] const std::array<double,
                                 kDimension> &GetRelativePosition() const;


  inline void writeConfig(std::ofstream &ofs) const;
  inline bool readPOSCAR(std::ifstream &ifs, const bool &realOption);
  inline void writePrl(std::ofstream &ofs) const;

 private:
  // atom id which is an unique int for every atom indexed form 0
  int id_;
  double mass_;
  std::string type_;
  // absolute position
  std::array<double, kDimension> absolute_position_{};
  // relative position in the box
  std::array<double, kDimension> relative_position_{};
  // First and second nearest neighbor list
  std::array<int, kNumNearNeighbors> near_neighbor_list_{};
  // First nearest neighbor list
  std::array<int, kNumFirstNearestNeighbors> first_nearest_neighbor_list_{};

};

#endif //KN_INCLUDE_ATOM_H_
