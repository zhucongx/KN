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
      std::sort(atom_array_.begin(), atom_array_.end(), [](const auto &lhs, const auto &rhs) {
        // sort to make sure the uniqueness of clusters
        return lhs.GetId() < rhs.GetId();
      });
    }
    template<typename ... Ts,
        typename std::enable_if<sizeof...(Ts) == DataSize, int>::type = 0>
    explicit Cluster(Ts &&... ts) : atom_array_{std::forward<Ts>(ts)...} {
      std::sort(atom_array_.begin(), atom_array_.end(), [](const auto &lhs, const auto &rhs) {
        // sort to make sure the uniqueness of clusters
        return lhs.GetId() < rhs.GetId();
      });
    }

    /// Getter
    [[nodiscard]] const Atom &GetAtomAt(size_t i) const {
      return atom_array_[i];
    }
    // [[nodiscard]] const Atom &GetIdAt(size_t i) const {
    //   return atom_array_[i].GetId() ;
    // }
    ///Operators
    friend bool operator==(const Cluster<DataSize> &lhs, const Cluster<DataSize> &rhs) {
      for (size_t i = 0; i < DataSize; ++i) {
        if (lhs.atom_array_[i].GetId() != rhs.atom_array_[i].GetId())
          return false;
      }
      return true;
    }
    friend bool operator<(const Cluster<DataSize> &lhs, const Cluster<DataSize> &rhs) {
      for (size_t i = 0; i < DataSize; ++i) {
        if (lhs.atom_array_[i].GetId() < rhs.atom_array_[i].GetId())
          return true;
      }
      return false;
    }
    friend size_t hash_value(const Cluster<DataSize> &cluster) {
      size_t seed = 0;
      for (size_t i = 0; i < DataSize; ++i) {
        boost::hash_combine(seed, cluster.atom_array_[i].GetId());
      }
      return seed;
    }
  private:
    std::array<Atom, DataSize> atom_array_;
};

}// namespace cfg


#endif //KN_KN_CFG_INCLUDE_CLUSTER_HPP_
