#ifndef KN_SRC_BONDCOUNTER_H_
#define KN_SRC_BONDCOUNTER_H_

#include <set>
#include "Config.h"

namespace box {

class BondCounter {
 public:

  BondCounter();
  BondCounter(Vector3<double> factor,
              Vector3<double> plane,
              Vector3<double> burger_vector);
  void SetFactor(const Vector3<double> &factor);
  void SetPlane(Vector3<double> miller_index);
  void SetBurgersVector(Vector3<double> miller_index);
  [[nodiscard]] std::map<Bond, int> GetBondChange() const;
  Config config_;
 private:
  // make sure d1 < d2
  static std::vector<Atom::Rank> GetAtomListBetweenPlanesHelper(
      const Config &config,
      const Vector3<double> &abc,
      const double &d1,
      const double &d2);
  static std::map<Bond, int> CountBondsBetweenTwoGroupHelper(
      const Config &config,
      const std::vector<Atom::Rank> &group1,
      const std::vector<Atom::Rank> &group2);
  [[nodiscard]] Vector3<double> GetPlaneDistanceVector(const Vector3<double> &plane_index) const;
  [[nodiscard]] Vector3<double> GetBurgerDistanceVector(const Vector3<double> &burger_vector) const;
  std::set<Vector3<double>> plane_set_;
  std::set<Vector3<double>> burgers_vector_set_;
  Vector3<double> factor_{};
};

}// namespace box
#endif //KN_SRC_BONDCOUNTER_H_
