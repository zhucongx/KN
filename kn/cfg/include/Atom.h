#ifndef KN_KN_CFG_INCLUDE_ATOM_H_
#define KN_KN_CFG_INCLUDE_ATOM_H_
#include <fstream>
#include <array>
#include <vector>
#include <unordered_set>

#include "VectorMatrix.hpp"
#include "Config.h"
namespace cfg {
class Config;
class Atom {
  public:
    /// Constructor
    Atom();
    // Set both relative and absolute position, but will be corrected later
    Atom(size_t id, double mass, std::string type, double x, double y, double z);
    Atom(size_t id, double mass, std::string type, Vector_t position);
    /// Getter
    [[nodiscard]] const Vector_t &GetCartesianPosition() const;
    [[nodiscard]] const Vector_t &GetRelativePosition() const;
    [[nodiscard]] size_t GetId() const;
    [[nodiscard]] double GetMass() const;
    [[nodiscard]] const std::string &GetType() const;
    [[nodiscard]] const std::vector<size_t> &GetFirstNearestNeighborsList() const;
    [[nodiscard]] const std::vector<size_t> &GetSecondNearestNeighborsList() const;
    [[nodiscard]] const std::vector<size_t> &GetThirdNearestNeighborsList() const;
    [[nodiscard]] std::unordered_set<size_t> GetFirstAndSecondThirdNeighborsSet() const;
    /// Operators
    // bool operator<(const Atom &rhs) const;
    /// Setter
    void SetCartesianPosition(const Vector_t &cartesian_position);
    void SetRelativePosition(const Vector_t &relative_position);
    void SetId(size_t id);
    void SetType(const std::string &type);
    /// Function related to neighbors list
    void AppendFirstNearestNeighborsList(size_t index);
    void AppendSecondNearestNeighborsList(size_t index);
    void AppendThirdNearestNeighborsList(size_t index);
    void CleanNeighborsLists();
  private:
    size_t id_{};
    double mass_{};
    std::string type_;
    // absolute position
    Vector_t cartesian_position_{};
    // relative position in the box
    Vector_t relative_position_{};
    // // near neighbor hashset
    // std::unordered_set<size_t> near_neighbor_hashset_;
    // First nearest neighbor list
    std::vector<size_t> first_nearest_neighbors_list_;
    // Second nearest neighbor list
    std::vector<size_t> second_nearest_neighbors_list_;
    // Third nearest neighbor list
    std::vector<size_t> third_nearest_neighbors_list_;
    // Fourth nearest neighbor list
    // std::vector<size_t> fourth_nearest_neighbors_list_;
    // atom id which is an unique Rank for every atom indexed form 0
  public:
    /// Friend function
    friend void AtomsJump(Config &config, size_t lhs, size_t rhs);
};
Vector_t GetRelativeDistanceVector(const Atom &first, const Atom &second);
double FindMass(const std::string &elem);

}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_ATOM_H_
