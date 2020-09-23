#ifndef KN_KN_CFG_INCLUDE_CLUSTERS_HPP_
#define KN_KN_CFG_INCLUDE_CLUSTERS_HPP_
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include "Atom.h"
namespace cfg {
class Atom;

template<size_t DataSize>
class Clusters {
  public:
    /// Constructor
    explicit Clusters(const std::array<Atom, DataSize> &atom_array);
    /// Getter
    [[nodiscard]] int GetIndexAt(int i) const;
    ///Operators
    friend bool operator==(const Clusters<DataSize> &lhs, const Clusters<DataSize> &rhs);
    friend bool operator<(const Clusters<DataSize> &lhs, const Clusters<DataSize> &rhs);
    friend size_t hash_value(const Clusters<DataSize> &cluster);
  private:
    std::array<int, DataSize> atom_array_;
};

template<size_t DataSize>
Clusters<DataSize>::Clusters(const std::array<Atom, DataSize> &atom_array)
    :atom_array_(atom_array) {
  std::sort(atom_array.begin(), atom_array.end(), [](int a, int b) {
    // sort to make sure the uniqueness of clusters
    return a < b;
  });
}
template<size_t DataSize>
int Clusters<DataSize>::GetIndexAt(int i) const{
  return atom_array_[i];
}

template<size_t DataSize>
bool operator==(const Clusters<DataSize> &lhs, const Clusters<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    if (lhs.atom_array_[i] != rhs.atom_array_[i])
      return false;
  }
  return true;
}
template<size_t DataSize>
bool operator<(const Clusters<DataSize> &lhs, const Clusters<DataSize> &rhs) {
  for (size_t i = 0; i < DataSize; ++i) {
    if (lhs.atom_array_[i] < rhs.atom_array_[i])
      return true;
  }
  return false;
}
template<size_t DataSize>
size_t hash_value(const Clusters<DataSize> &cluster) {
  size_t seed = 0;
  for (size_t i = 0; i < DataSize; ++i) {
    boost::hash_combine(seed, cluster.atom_array_[i]);
  }
  return seed;
}
}// namespace cfg


#endif //KN_KN_CFG_INCLUDE_CLUSTERS_HPP_
