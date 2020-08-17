#ifndef KN_SRC_BOND_H_
#define KN_SRC_BOND_H_
#include <iostream>

namespace kn {
class Bond {
  public:
    Bond(std::string type1, std::string type2);
    [[nodiscard]] const std::string &GetType1() const;
    [[nodiscard]] const std::string &GetType2() const;
    friend std::ostream &operator<<(std::ostream &os, const Bond &bond);
    friend bool operator<(const Bond &lhs, const Bond &rhs);
    friend bool operator==(const Bond &lhs, const Bond &rhs);
    friend bool operator!=(const Bond &lhs, const Bond &rhs);
  private:
    std::string type1_;
    std::string type2_;
};
} // namespace kn

namespace std {
template<>
struct hash<kn::Bond> {
  std::size_t operator()(const kn::Bond &k) const {
    return (std::hash<std::string>()(k.GetType1())) ^ (std::hash<std::string>()(k.GetType2()));
  }
};
}

#endif //KN_SRC_BOND_H_
