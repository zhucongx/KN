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
    Atom();
    // Set both relative and absolute position, but will be corrected later
    Atom(int id, double mass, std::string type, double x, double y, double z);
    Atom(int id, double mass, std::string type, Vector3 position);

  public:
    // absolute position
    Vector3 cartesian_position_{};
    // relative position in the box
    Vector3 relative_position_{};
    // Second nearest neighbor list
    std::vector<int> second_nearest_neighbor_list_;
    // First nearest neighbor list
    std::vector<int> first_nearest_neighbor_list_;
    // atom id which is an unique Rank for every atom indexed form 0
    int id_{};
    double mass_{};
    std::string type_;
};
} // namespace kn
#endif //KN_SRC_ATOM_H_
