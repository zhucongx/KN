#ifndef KN_KN_CFG_INCLUDE_BOND_H_
#define KN_KN_CFG_INCLUDE_BOND_H_
#include <iostream>
#include <boost/functional/hash.hpp>

namespace cfg {
class Bond {
  public:
    /// Constructor
    Bond(size_t label, std::string type1, std::string type2)
        : label_(label), type1_(std::move(type1)), type2_(std::move(type2)) {
      if (type2_.compare(type1_) < 0)
        std::swap(type1_, type2_);
    }
    /// Getter
    [[nodiscard]] const std::string &GetType1() const {
      return type1_;
    }
    [[nodiscard]] const std::string &GetType2() const {
      return type2_;
    }
    /// Operators
    friend std::ostream &operator<<(std::ostream &os, const Bond &bond) {
      os << bond.label_ << bond.type1_ << "-" << bond.type2_;
      return os;
    }
    friend bool operator<(const Bond &lhs, const Bond &rhs) {
      if (lhs.label_ < rhs.label_)
        return true;
      if (rhs.label_ < lhs.label_)
        return false;
      if (lhs.type1_ < rhs.type1_)
        return true;
      if (rhs.type1_ < lhs.type1_)
        return false;
      return lhs.type2_ < rhs.type2_;
    }
    friend bool operator==(const Bond &lhs, const Bond &rhs) {
      return lhs.label_ == rhs.label_ && lhs.type1_ == rhs.type1_ && lhs.type2_ == rhs.type2_;
    }
    friend size_t hash_value(const Bond &bond) {
      size_t seed = 0;
      boost::hash_combine(seed, bond.label_);
      boost::hash_combine(seed, bond.type1_);
      boost::hash_combine(seed, bond.type2_);
      return seed;
    }
  private:
    size_t label_;
    std::string type1_;
    std::string type2_;
};
}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_BOND_H_
