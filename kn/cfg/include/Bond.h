#ifndef KN_KN_CFG_INCLUDE_BOND_H_
#define KN_KN_CFG_INCLUDE_BOND_H_
#include <iostream>

#include <boost/functional/hash.hpp>

namespace cfg {
class Bond {
  public:
    /// Constructor
    Bond(std::string type1, std::string type2);
    /// Getter
    [[nodiscard]] const std::string &GetType1() const;
    [[nodiscard]] const std::string &GetType2() const;
    /// Operators
    friend std::ostream &operator<<(std::ostream &os, const Bond &bond);
    friend bool operator<(const Bond &lhs, const Bond &rhs);
    friend bool operator==(const Bond &lhs, const Bond &rhs);
    friend size_t hash_value(const Bond &bond);
  private:
    std::string type1_;
    std::string type2_;
};
}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_BOND_H_
