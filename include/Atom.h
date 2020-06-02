#ifndef KN_SRC_ATOM_H_
#define KN_SRC_ATOM_H_

#include <fstream>
#include <array>
#include <vector>

#include "Constants.h"
#include "VectorMatrix.h"

namespace kn {

class Atom {
 public:
  typedef int Rank;

  Atom();
  // Set both relative and absolute position, but will be corrected later
  Atom(Rank id, double mass, std::string type, double x, double y, double z);
  Atom(Rank id, double mass, std::string type, Vector3 position);

 public:
  // absolute position
  Vector3 cartesian_position_{};
  // relative position in the box
  Vector3 relative_position_{};
  // First and second nearest neighbor list
  std::vector<Rank> second_nearest_neighbor_list_;
  // First nearest neighbor list
  std::vector<Rank> first_nearest_neighbor_list_;
  // atom id which is an unique Rank for every atom indexed form 0
  Rank id_{};
  double mass_{};
  std::string type_;
};

Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second);
}// namespace kn
#endif //KN_SRC_ATOM_H_
