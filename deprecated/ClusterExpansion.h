#ifndef KN_KN_KMC_INCLUDE_CLUSTEREXPANSION_H_
#define KN_KN_KMC_INCLUDE_CLUSTEREXPANSION_H_

#include <utility>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "Config.h"

namespace kn::ClusterExpansion {
// Todo rewrite this namespace
// initial and final state
// this function rotate and sort the config in a particular way, basically from center to outside
std::vector<double> GetAverageClusterFunctions(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

std::vector<double> GetAverageClusterFunctionsBack(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

class Pair {
  public:
    Pair(cfg::Atom atom1, cfg::Atom atom2) : atom1_(std::move(atom1)), atom2_(std::move(atom2)) {
      if (atom2_.GetId() < atom1_.GetId()) std::swap(atom1_, atom2_);
    }
    bool operator==(const Pair &rhs) const {
      return atom1_.GetId() == rhs.atom1_.GetId() &&
          atom2_.GetId() == rhs.atom2_.GetId();
    }
    bool operator<(const Pair &rhs) const {
      if (atom1_.GetId() < rhs.atom1_.GetId())
        return true;
      if (rhs.atom1_.GetId() < atom1_.GetId())
        return false;
      return atom2_.GetId() < rhs.atom2_.GetId();
    }
    [[nodiscard]] const cfg::Atom &GetAtom1() const {
      return atom1_;
    }
    [[nodiscard]] const cfg::Atom &GetAtom2() const {
      return atom2_;
    }
    friend size_t hash_value(Pair const &pair) {
      size_t seed = 0;
      boost::hash_combine(seed, pair.atom1_.GetId());
      boost::hash_combine(seed, pair.atom2_.GetId());
      return seed;
    }

  private:
    cfg::Atom atom1_;
    cfg::Atom atom2_;
};

class Triplet {
  public:
    Triplet(cfg::Atom atom1, cfg::Atom atom2, cfg::Atom atom3)
        : atom1_(std::move(atom1)), atom2_(std::move(atom2)), atom3_(std::move(atom3)) {
      if (atom2_.GetId() < atom1_.GetId()) std::swap(atom1_, atom2_);
      if (atom3_.GetId() < atom2_.GetId()) std::swap(atom2_, atom3_);
      if (atom2_.GetId() < atom1_.GetId()) std::swap(atom1_, atom2_);
    }
    bool operator==(const Triplet &rhs) const {
      return atom1_.GetId() == rhs.atom1_.GetId() &&
          atom2_.GetId() == rhs.atom2_.GetId() &&
          atom3_.GetId() == rhs.atom3_.GetId();
    }
    bool operator<(const Triplet &rhs) const {
      if (atom1_.GetId() < rhs.atom1_.GetId())
        return true;
      if (rhs.atom1_.GetId() < atom1_.GetId())
        return false;
      if (atom2_.GetId() < rhs.atom2_.GetId())
        return true;
      if (rhs.atom2_.GetId() < atom2_.GetId())
        return false;
      return atom3_.GetId() < rhs.atom3_.GetId();
    }
    [[nodiscard]] const cfg::Atom &GetAtom1() const {
      return atom1_;
    }
    [[nodiscard]] const cfg::Atom &GetAtom2() const {
      return atom2_;
    }
    [[nodiscard]] const cfg::Atom &GetAtom3() const {
      return atom3_;
    }
    friend size_t hash_value(Triplet const &triplet) {
      size_t seed = 0;
      boost::hash_combine(seed, triplet.atom1_.GetId());
      boost::hash_combine(seed, triplet.atom2_.GetId());
      boost::hash_combine(seed, triplet.atom3_.GetId());
      return seed;
    }

  private:
    cfg::Atom atom1_;
    cfg::Atom atom2_;
    cfg::Atom atom3_;
};

} // namespace kn::ClusterExpansion
#endif //KN_KN_KMC_INCLUDE_CLUSTEREXPANSION_H_
