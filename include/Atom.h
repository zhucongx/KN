#ifndef KN_SRC_ATOM_H_
#define KN_SRC_ATOM_H_

#include <fstream>
#include <array>
#include <vector>

#include "Constants.h"
#include "Matrix33.h"

namespace box
{

class Atom
{
 public:
  typedef int Rank;

  Atom();
  // Set both relative and absolute position, but will be corrected later
  Atom(Rank id, double mass, std::string type, double x, double y, double z);
  Atom(Rank id, double mass, std::string type, Vector3<double> position);
  void SetId(Rank id);
  [[nodiscard]] Rank GetId() const;
  void SetMass(double mass);
  [[nodiscard]] double GetMass() const;
  void SetType(const std::string &type);
  [[nodiscard]] const std::string &GetType() const;

  // absolute position
  Vector3<double> cartesian_position_{};
  // relative position in the box
  Vector3<double> relative_position_{};
  // First and second nearest neighbor list
  std::vector<Rank> second_nearest_neighbor_list_;
  // First nearest neighbor list
  std::vector<Rank> first_nearest_neighbor_list_;
 private:
  // atom id which is an unique Rank for every atom indexed form 0
  Rank id_{};
  double mass_{};
  std::string type_;
};

Vector3<double> GetRelativeDistanceVector(const Atom &first, const Atom &second);
}// namespace box
#endif //KN_SRC_ATOM_H_
