#ifndef KN_SRC_BOND_H_
#define KN_SRC_BOND_H_
#include <iostream>
class Bond {
 public:
  Bond(std::string type1, std::string type_2);
  bool operator<(const Bond &rhs) const;
  friend std::ostream &operator<<(std::ostream &os, const Bond &bond);
 private:
  std::string type1_;
  std::string type2_;
};

#endif //KN_SRC_BOND_H_
