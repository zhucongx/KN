#include "BondCounter.h"
namespace box {
BondCounter::BondCounter() = default;
BondCounter::BondCounter(Vector3<double> factor,
                         Vector3<double> plane,
                         Vector3<double> burger_vector) : factor_(factor) {
  SetPlane(plane);
  SetBurgersVector(burger_vector);
}
void BondCounter::SetFactor(const Vector3<double> &factor) {
  factor_ = factor;
}
void BondCounter::SetPlane(Vector3<double> miller_index) {
  for (const double &i :{-1.0, 1.0}) {
    for (const double &j :{-1.0, 1.0}) {
      for (const double &k :{-1.0, 1.0}) {
        for (int l = 0; l < 3; l++) {
          miller_index = {miller_index.y, miller_index.z, miller_index.x};
          plane_set_.insert(StarProduct(miller_index, {i, j, k}));
        }
      }
    }
  }
}
void BondCounter::SetBurgersVector(Vector3<double> miller_index) {
  for (int l = 0; l < 3; l++) {
    miller_index = {miller_index.y, miller_index.z, miller_index.x};
    for (const double &i :{-1.0, 1.0}) {
      for (const double &j :{-1.0, 1.0}) {
        for (const double &k :{-1.0, 1.0}) {
          burgers_vector_set_.insert(StarProduct(miller_index, {i, j, k}));
        }
      }
    }
  }
}

std::map<Bond, int> BondCounter::GetBondChange() const {
  std::map<Bond, int> bonds_changed;
  for (const auto &plane_index:plane_set_) {
    auto unsliped_config = config_;
    Vector3<double>
        plane_distance_vector = GetPlaneDistanceVector(plane_index);
    /// {1/90, 1/90, 1/90}
    unsliped_config.MoveRelativeDistance(0.5 * plane_distance_vector);
    // How many space in a unit cell does a plane divide
    // {100}1 {110}2 {111}3 {200}4 {220}4 {222}6
    // d_spacing is d=1/sqrt(h^2+k^2+l^2)
    // the length along that direction is L=h'^2+k'^2+l'2
    // where (h',k',l') is (h,k,l)/GCD(hkl)
    // if n is the number of space divided in one cubic
    // n = L/d_spacing = (h^2+k^2+l^2)/GCD(hkl)
    // so because of periodic boundary condition. We only need to iterate
    // LCM(1/n*factor) times
    /// 30 here for FCC
    int iteration_time = 30;

    const double d1 = 1.5;
    const double d2 = d1 - Sum(plane_distance_vector);
    const double d3 = d1 + Sum(plane_distance_vector);
    for (int i = 0; i < iteration_time; i++) {
      unsliped_config.MoveRelativeDistance(plane_distance_vector);
      // We select four planes, and we move two plan by the burgers vector and
      // calculate the bonds change between the these two planes and the
      // other two.
      std::vector<Atom::Rank>
          atoms_on_plane1 = GetAtomListBetweenPlanesHelper(unsliped_config,
                                                           plane_index,
                                                           d2, d1);
      std::vector<Atom::Rank>
          atoms_on_plane2 = GetAtomListBetweenPlanesHelper(unsliped_config,
                                                           plane_index,
                                                           d1, d3);

      std::map<Bond, int>
          bonds_map_before = CountBondsBetweenTwoGroupHelper(unsliped_config,
                                                             atoms_on_plane1,
                                                             atoms_on_plane2);
      //  move atoms
      for (const auto &burgers_vector:burgers_vector_set_) {
        if (DotProduct(plane_index, burgers_vector) != 0)
          continue;
        auto burger_distance_vector = GetBurgerDistanceVector(burgers_vector);
        Config sliped_config = unsliped_config;
        for (const auto &index:atoms_on_plane2) {
          sliped_config
              .MoveOneAtomRelativeDistance(index, burger_distance_vector);
        }
        std::map<Bond, int>
            bonds_map_after = CountBondsBetweenTwoGroupHelper(sliped_config,
                                                              atoms_on_plane1,
                                                              atoms_on_plane2);
        // bonds_map_temp = bonds_map_after - bonds_map_before;
        for (const auto&[key, count] : bonds_map_after) {
          bonds_changed[key] += count;
        }
        for (const auto&[key, count] : bonds_map_before) {
          bonds_changed[key] -= count;
        }
      }
    }
  }
  return bonds_changed;
}

// relative distance
Vector3<double> BondCounter::GetPlaneDistanceVector(const Vector3<double> &plane_index) const {
  double inner = InnerProduct(plane_index);
  Vector3<double> distance_index = (1.0 / inner * plane_index);
  return StarDivide(distance_index, factor_);
}

Vector3<double> BondCounter::GetBurgerDistanceVector(const Vector3<double> &burger_vector) const {
  return StarDivide(burger_vector, factor_);
}

std::vector<Atom::Rank> BondCounter::GetAtomListBetweenPlanesHelper(
    const Config &config,
    const Vector3<double> &abc,
    const double &d1,
    const double &d2) {
  std::vector<Atom::Rank> rank_list;
  double larger = std::max(d1, d2);
  double smaller = std::min(d1, d2);
  for (Atom::Rank i = 0; i < config.GetNumAtoms(); ++i) {
    double d_checked = DotProduct(config.GetAtom(i).relative_position_, abc);
    if (d_checked > larger || d_checked < smaller)
      continue;
    rank_list.push_back(i);
  }
  return rank_list;
}
std::map<Bond, int> BondCounter::CountBondsBetweenTwoGroupHelper(
    const Config &config,
    const std::vector<Atom::Rank> &group1,
    const std::vector<Atom::Rank> &group2) {
  std::map<Bond, int> map_out;
  for (const auto &index1:group1) {
    for (const auto &index2:group2) {
      Vector3<double> relative_distance_vector =
          GetRelativeDistanceVector(config.GetAtom(index1),
                                    config.GetAtom(index2));
      double absolute_distance_squared = (InnerProduct(
          relative_distance_vector * config.GetBravaisMatrix()));

      if (absolute_distance_squared < 9) {
        map_out[{config.GetAtom(index1).GetType(),
                 config.GetAtom(index2).GetType()}]++;
      }
    }
  }
  return map_out;
}

}// namespace box