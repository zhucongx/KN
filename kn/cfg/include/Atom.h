#ifndef KN_SRC_ATOM_H_
#define KN_SRC_ATOM_H_

#include <fstream>
#include <array>
#include <vector>
#include <unordered_set>

#include "Constants.h"
#include "VectorMatrix.h"

namespace cfg {
class Atom;
Vector3 GetRelativeDistanceVector(const Atom &first, const Atom &second);
void AtomsJump(Atom &lhs, Atom &rhs);


class Atom {
  public:
    Atom();
    // Set both relative and absolute position, but will be corrected later
    Atom(int id, double mass, std::string type, double x, double y, double z);
    Atom(int id, double mass, std::string type, Vector3 position);
    friend bool operator==(const Atom &lhs, const Atom &rhs);
    friend bool operator<(const Atom &lhs, const Atom &rhs);
    [[nodiscard]] const Vector3 &GetCartesianPosition() const;
    void SetCartesianPosition(const Vector3 &cartesian_position);
    [[nodiscard]] const Vector3 &GetRelativePosition() const;
    void SetRelativePosition(const Vector3 &relative_position);
    [[nodiscard]] int GetId() const;
    void SetId(int id);
    [[nodiscard]] double GetMass() const;
    [[nodiscard]] const std::string &GetType() const;
    void SetType(const std::string &type);
    [[nodiscard]] const std::vector<int> &GetFirstNearestNeighborsList() const;
    [[nodiscard]] const std::vector<int> &GetSecondNearestNeighborsList() const;
    [[nodiscard]] const std::vector<int> &GetThirdNearestNeighborsList() const;
    [[nodiscard]]  std::unordered_set<int> GetFirstAndSecondNeighborsSet() const;
    [[nodiscard]]  std::unordered_set<int> GetFirstAndSecondThirdNeighborsSet() const;

    void AppendFirstNearestNeighborsList(int index);
    void AppendSecondNearestNeighborsList(int index);
    void AppendThirdNearestNeighborsList(int index);

    void CleanNeighborsLists();

    friend void cfg::AtomsJump(Atom &lhs, Atom &rhs);
  private:
    int id_{};
    double mass_{};
    std::string type_;

    // absolute position
    Vector3 cartesian_position_{};
    // relative position in the box
    Vector3 relative_position_{};
    // near neighbor hashset
    std::unordered_set<int> near_neighbor_hashset_;
    // First nearest neighbor list
    std::vector<int> first_nearest_neighbors_list_;
    // Second nearest neighbor list
    std::vector<int> second_nearest_neighbors_list_;
    // Third nearest neighbor list
    std::vector<int> third_nearest_neighbors_list_;
    // Fourth nearest neighbor list
    // std::vector<int> fourth_nearest_neighbors_list_;
    // atom id which is an unique Rank for every atom indexed form 0
};
} // namespace cfg
#endif //KN_SRC_ATOM_H_
