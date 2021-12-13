#ifndef KN_KN_CFG_INCLUDE_CLUSTER_HPP_
#define KN_KN_CFG_INCLUDE_CLUSTER_HPP_
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include<type_traits>
#include <boost/functional/hash.hpp>
#include "Atom.h"
namespace cfg {
class Atom;

template<size_t DataSize>
class Cluster {
  public:
    /// Constructor
    explicit Cluster(std::array<Atom, DataSize> atom_array) : atom_array_(std::move(atom_array)) {
      Sort();
    }
    template<typename ... Ts, typename std::enable_if<sizeof...(Ts) == DataSize, int>::type = 0>
    explicit Cluster(Ts &&... ts) : atom_array_{std::forward<Ts>(ts)...} {
      Sort();
    }

    /// Getter
    [[nodiscard]] const Atom &GetAtomAt(size_t i) const {
      return atom_array_[i];
    }

    ///Operators
    friend bool operator==(const Cluster<DataSize> &lhs, const Cluster<DataSize> &rhs) {
      for (size_t i = 0; i < DataSize; ++i) {
        if (lhs.atom_array_[i].GetId() != rhs.atom_array_[i].GetId())
          return false;
      }
      return true;
    }
    friend size_t hash_value(const Cluster<DataSize> &cluster) {
      size_t seed = 0;
      for (size_t i = 0; i < DataSize; ++i) {
        boost::hash_combine(seed, cluster.atom_array_[i].GetId());
      }
      return seed;
    }
  private:
    void Sort() {
      std::sort(atom_array_.begin(), atom_array_.end(), [](const auto &lhs, const auto &rhs) {
        const auto &relative_position_lhs = lhs.GetRelativePosition();
        const auto &relative_position_rhs = rhs.GetRelativePosition();
        const double diff_norm =
            Inner(relative_position_lhs - 0.5) - Inner(relative_position_rhs - 0.5);
        if (diff_norm < -kEpsilon) { return true; }
        if (diff_norm > kEpsilon) { return false; }
        const double diff_x = std::abs(relative_position_lhs[kXDimension] - 0.5)
            - std::abs(relative_position_rhs[kXDimension] - 0.5);
        if (diff_x < -kEpsilon) { return true; }
        if (diff_x > kEpsilon) { return false; }
        const double diff_x2 =
            relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
        if (diff_x2 < -kEpsilon) { return true; }
        if (diff_x2 > kEpsilon) { return false; }
        const double diff_y =
            relative_position_lhs[kYDimension] - relative_position_rhs[kYDimension];
        if (diff_y < -kEpsilon) { return true; }
        if (diff_y > kEpsilon) { return false; }
        const double diff_z =
            relative_position_lhs[kZDimension] - relative_position_rhs[kZDimension];
        if (diff_z < -kEpsilon) { return true; }
        if (diff_z > kEpsilon) { return false; }
        return lhs.GetId() < rhs.GetId();
      });
    }
    std::array<Atom, DataSize> atom_array_;
};

}// namespace cfg


#endif //KN_KN_CFG_INCLUDE_CLUSTER_HPP_
