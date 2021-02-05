#include "Bond.h"
namespace cfg {
Bond::Bond(std::string type1, std::string type2)
    : type1_(std::move(type1)), type2_(std::move(type2)) {
  if (type2_.compare(type1_) < 0)
    std::swap(type1_, type2_);
}
const std::string &Bond::GetType1() const {
  return type1_;
}
const std::string &Bond::GetType2() const {
  return type2_;
}
std::ostream &operator<<(std::ostream &os, const Bond &bond) {
  os << bond.type1_ << "-" << bond.type2_;
  return os;
}
bool operator<(const Bond &lhs, const Bond &rhs) {
  if (lhs.type1_ < rhs.type1_)
    return true;
  if (rhs.type1_ > lhs.type1_)
    return false;
  return lhs.type2_ < rhs.type2_;
}
bool operator==(const Bond &lhs, const Bond &rhs) {
  return lhs.type1_ == rhs.type1_ &&
      lhs.type2_ == rhs.type2_;
}
size_t hash_value(const Bond &bond) {
  size_t seed = 0;
  boost::hash_combine(seed, bond.type1_);
  boost::hash_combine(seed, bond.type2_);
  return seed;
}
}// namespace cfg