#ifndef KN_INCLUDE_BONDCOUNTER_H_
#define KN_INCLUDE_BONDCOUNTER_H_

#include <set>
#include "Config.h"

namespace box {

class BondCounter {
 public:
  BondCounter(const std::string& cfg_file_name,
              Vector3 factor,
              Vector3 plane,
              Vector3 burger_vector);
  void SetFactor(const Vector3 &factor);
  void SetPlane(Vector3 miller_index);
  void SetBurgersVector(Vector3 miller_index);
  [[nodiscard]] std::map<Bond, int> GetBondChange() const;
 private:
  static double FindD3Helper(const Vector3 &plane_index,
                             const Vector3 &box_low_bound,
                             const Vector3 &box_high_bound);
  static std::map<Bond, int> CountBondsBetweenTwoGroupHelper(
      const Config &config,
      const std::vector<Atom::Rank> &group1,
      const std::vector<Atom::Rank> &group2);
  [[nodiscard]] Vector3 GetPlaneDistanceVector(const Vector3 &plane_index) const;
  [[nodiscard]] Vector3 GetBurgerDistanceVector(const Vector3 &burger_vector) const;

  Config config_;
  std::set<Vector3> plane_set_;
  std::set<Vector3> burgers_vector_set_;
  Vector3 factor_{};
};

}// namespace box
#endif //KN_INCLUDE_BONDCOUNTER_H_
