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
    Atom(int id, double mass, std::string type, double x, double y, double z);
    Atom(int id, double mass, std::string type, Vector_t position);
    /// Getter
    [[nodiscard]] const Vector_t &GetCartesianPosition() const;
    [[nodiscard]] const Vector_t &GetRelativePosition() const;
    [[nodiscard]] int GetId() const;
    [[nodiscard]] double GetMass() const;
    [[nodiscard]] const std::string &GetType() const;
    [[nodiscard]] const std::vector<int> &GetFirstNearestNeighborsList() const;
    [[nodiscard]] const std::vector<int> &GetSecondNearestNeighborsList() const;
    [[nodiscard]] const std::vector<int> &GetThirdNearestNeighborsList() const;
    [[nodiscard]] std::unordered_set<int> GetFirstAndSecondThirdNeighborsSet() const;
    /// Operators
    // bool operator<(const Atom &rhs) const;
    /// Setter
    void SetCartesianPosition(const Vector_t &cartesian_position);
    void SetRelativePosition(const Vector_t &relative_position);
    void SetId(int id);
    void SetType(const std::string &type);
    /// Function related to neighbors list
    void AppendFirstNearestNeighborsList(int index);
    void AppendSecondNearestNeighborsList(int index);
    void AppendThirdNearestNeighborsList(int index);
    void CleanNeighborsLists();
  private:
    int id_{};
    double mass_{};
    std::string type_;
    // absolute position
    Vector_t cartesian_position_{};
    // relative position in the box
    Vector_t relative_position_{};
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
  public:
    /// Friend function
    friend void AtomsJump(Config &config, int lhs, int rhs);
};
Vector_t GetRelativeDistanceVector(const Atom &first, const Atom &second);
double FindMass(const std::string &elem);

}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_ATOM_H_
