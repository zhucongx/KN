#include "Bond.h"
namespace kn {

Bond::Bond(std::string type1, std::string type2)
  : type1_(std::move(type1)), type2_(std::move(type2)) {
  if (type2_.compare(type1_) < 0)
    std::swap(type1_, type2_);
}

bool Bond::operator<(const Bond &rhs) const {
  if (type1_ < rhs.type1_)
    return true;
  if (rhs.type1_ < type1_)
    return false;
  return type2_ < rhs.type2_;
}

std::ostream &operator<<(std::ostream &os, const Bond &bond) {
  os << bond.type1_ << "-" << bond.type2_;
  return os;
}
} // namespace kn
