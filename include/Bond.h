#ifndef KN_SRC_BOND_H_
#define KN_SRC_BOND_H_
#include <iostream>
namespace kn {
class Bond {
 public:
  Bond(std::string type1, std::string type2);
  bool operator<(const Bond &rhs) const;
  friend std::ostream &operator<<(std::ostream &os, const Bond &bond);
 private:
  std::string type1_;
  std::string type2_;
};
}// namespace kn
#endif //KN_SRC_BOND_H_
