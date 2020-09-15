#ifndef KN_SRC_BOND_H_
#define KN_SRC_BOND_H_
#include <iostream>

#include <boost/functional/hash.hpp>

namespace cfg {



class Bond {
  public:
    Bond(std::string type1, std::string type2)
        : type1_(std::move(type1)), type2_(std::move(type2)) {
      if (type2_.compare(type1_) < 0)
        std::swap(type1_, type2_);
    }
    [[nodiscard]] const std::string &GetType1() const {
      return type1_;
    }
    [[nodiscard]] const std::string &GetType2() const {
      return type2_;
    }
    friend std::ostream &operator<<(std::ostream &os, const Bond &bond) {
      os << bond.type1_ << "-" << bond.type2_;
      return os;
    }
    friend bool operator<(const Bond &lhs, const Bond &rhs) {
      if (lhs.type1_ < rhs.type1_)
        return true;
      if (rhs.type1_ < lhs.type1_)
        return false;
      return lhs.type2_ < rhs.type2_;
    }
    friend bool operator==(const Bond &lhs, const Bond &rhs) {
      return lhs.type1_ == rhs.type1_ &&
          lhs.type2_ == rhs.type2_;
    }
    friend std::size_t hash_value(Bond const &bond){
      std::size_t seed = 0;
      boost::hash_combine(seed,bond.type1_);
      boost::hash_combine(seed,bond.type2_);
            return seed;
    }
  private:
    std::string type1_;
    std::string type2_;
};



} // namespace cfg


// namespace std {
// template<>
// struct hash<cfg::Bond> {
//   // std::size_t operator()(const cfg::Bond &k) const {
//   //   return (std::hash<std::string>()(k.GetType1())) ^ (std::hash<std::string>()(k.GetType2()));
//   // }
//
//   std::size_t operator()(const cfg::Bond &k) const {
//     boost::hash<std::pair<std::string,std::string>> hasher;
//     return hasher({k.GetType1(),k.GetType2()});
//   }
// };
// }

#endif //KN_SRC_BOND_H_
