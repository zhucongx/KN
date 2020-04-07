#ifndef KN_SRC_ATOM_H_
#define KN_SRC_ATOM_H_

#include <fstream>
#include <array>

#include "ElemInfo.h"

// typedef int Rank;
// For FCC the first nearest neighbor is 12
constexpr int kNumFirstNearestNeighbors = 12;
// For FCC the first and second nearest neighbor is 18
constexpr int kNumNearNeighbors = 18;
// 3D simulation
constexpr int kDimension = 3;
enum : char { kXDim = 0, kYDim = 1, kZDim = 2 };

class Atom {
 public:
  explicit Atom(int id);
  Atom(int id, double mass, std::string type);
  Atom(int id, double mass, std::string type, double x, double y, double z);
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
#endif //KN_SRC_ATOM_H_
