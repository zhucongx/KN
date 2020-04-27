#ifndef KN_INCLUDE_BONDCOUNTER_H_
#define KN_INCLUDE_BONDCOUNTER_H_

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
 public:
  static double FindD3Helper(const Vector3<double> &plane_index,
                      const Vector3<double> &box_low_bound,
                      const Vector3<double> &box_high_bound);
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
#endif //KN_INCLUDE_BONDCOUNTER_H_
