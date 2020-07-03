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
    [[nodiscard]] const Vector3 &GetCartesianPosition() const;
    void SetCartesianPosition(const Vector3 &cartesian_position);
    [[nodiscard]] const Vector3 &GetRelativePosition() const;
    void SetRelativePosition(const Vector3 &relative_position);
    [[nodiscard]] const std::vector<int> &GetNearNeighborList() const;
    [[nodiscard]] const std::vector<int> &GetFirstNearestNeighborList() const;
    [[nodiscard]] int GetId() const;
    [[nodiscard]] double GetMass() const;
    [[nodiscard]] const std::string &GetType() const;
    void SetType(const std::string &type);
    void AppendNearNeighborList(int index);
    void AppendFirstNearestNeighborList(int index);

    friend Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second);
  private:
    // absolute position
    Vector3 cartesian_position_{};
    // relative position in the box
    Vector3 relative_position_{};
    // near neighbor list
    std::vector<int> near_neighbor_list_;
    // First nearest neighbor list
    std::vector<int> first_nearest_neighbor_list_;
    // atom id which is an unique Rank for every atom indexed form 0
    int id_{};
    double mass_{};
    std::string type_;
};
} // namespace kn
#endif //KN_SRC_ATOM_H_
