#ifndef KN_INCLUDE_CLUSTEREXPANSION_H_
#define KN_INCLUDE_CLUSTEREXPANSION_H_

#include <utility>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "Bond.h"
#include "Config.h"

namespace kn::ClusterExpansion {
// initial and final state
// this function rotate and sort the config in a particular way, basically from center to outside
std::vector<double> GetAverageClusterFunctions(
    const cfg::Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

std::vector<double> GetAverageClusterFunctionsBack(
    const cfg::Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, double> &type_category_hashmap);

class Pair {
  public:
    Pair(cfg::Atom atom1, cfg::Atom atom2) : atom1_(std::move(atom1)), atom2_(std::move(atom2)) {
      if (atom2_ < atom1_) std::swap(atom1_, atom2_);
    }
    bool operator==(const Pair &rhs) const {
      return atom1_ == rhs.atom1_ &&
          atom2_ == rhs.atom2_;
    }
    bool operator<(const Pair &rhs) const {
      if (atom1_ < rhs.atom1_)
        return true;
      if (rhs.atom1_ < atom1_)
        return false;
      return atom2_ < rhs.atom2_;
    }
    [[nodiscard]] const cfg::Atom &GetAtom1() const {
      return atom1_;
    }
    [[nodiscard]] const cfg::Atom &GetAtom2() const {
      return atom2_;
    }
    friend std::size_t hash_value(Pair const &pair) {
      std::size_t seed = 0;
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
      if (atom2_ < atom1_) std::swap(atom1_, atom2_);
      if (atom3_ < atom2_) std::swap(atom2_, atom3_);
      if (atom2_ < atom1_) std::swap(atom1_, atom2_);
    }
    bool operator==(const Triplet &rhs) const {
      return atom1_ == rhs.atom1_ &&
          atom2_ == rhs.atom2_ &&
          atom3_ == rhs.atom3_;
    }
    bool operator<(const Triplet &rhs) const {
      if (atom1_ < rhs.atom1_)
        return true;
      if (rhs.atom1_ < atom1_)
        return false;
      if (atom2_ < rhs.atom2_)
        return true;
      if (rhs.atom2_ < atom2_)
        return false;
      return atom3_ < rhs.atom3_;
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
    friend std::size_t hash_value(Triplet const &triplet) {
      std::size_t seed = 0;
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

}
#endif //KN_INCLUDE_CLUSTEREXPANSION_H_
