#ifndef KN_SRC_ATOM_H_
#define KN_SRC_ATOM_H_

#include <fstream>
#include <array>
#include <vector>

#include "Constants.h"
#include "VectorMatrix.h"

namespace kn {
struct Atom {
    Atom() = default;
    // Set both relative and absolute position, but will be corrected later
    Atom(int id, double mass, std::string type, double x, double y, double z) :
        cartesian_position_{x, y, z},
        relative_position_{x, y, z},
        id_(id),
        mass_(mass),
        type_(std::move(type)) {
    }

    Atom(int id, double mass, std::string type, Vector3 position) :
        cartesian_position_(position),
        relative_position_(position),
        id_(id),
        mass_(mass),
        type_(std::move(type)) {
    }

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
