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
  /// {1/90, 1/90, 1/90}
  // How many space in a unit cell does a plane divide
  // {100}1 {110}2 {111}3 {200}4 {220}4 {222}6
  // d_spacing is d=1/sqrt(h^2+k^2+l^2)
  // the length along that direction is L=h'^2+k'^2+l'2
  // where (h',k',l') is (h,k,l)/GCD(hkl)
  // if n is the number of space divided in one cubic
  // n = L/d_spacing = (h^2+k^2+l^2)/GCD(hkl)
  // select n planes if it's parallel to axes amd n-1 is not
  // so because of periodic boundary condition. We only need to iterate
  // LCM(1/n*factor) times
  /// 30 here for FCC
  int iteration_time = 30;
  /// slip plane should be d1 = 1/3(length_x+length_y+length_z) = 1, d2 = 2
  /// d1 ± 1/2Sum(plane_distance_vector) should help select one plane
  std::map<Bond, int> bonds_changed;
  for (const auto &plane_index:plane_set_) {
    auto unsliped_config = config_;
    Vector3<double> plane_distance_vector = GetPlaneDistanceVector(plane_index);

    const double delta = 0.5 * Sum(Abs(plane_distance_vector));
    double d3 = FindD3Helper(plane_index, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});

    const double d3_upper = d3 + delta;
    const double d3_lower = d3 - delta;
    const double d3_unslip_upper = d3 - delta;
    const double d3_unslip_lower = d3 - 3 * delta;

    const double d2 = d3 - 1;
    const double d2_upper = d2 + delta;
    const double d2_lower = d2 - delta;
    const double d2_unslip_upper = d2 - delta;
    const double d2_unslip_lower = d2 - 3 * delta;

    const double d1 = d3 - 2;
    const double d1_upper = d1 + delta;
    const double d1_lower = d1 - delta;
    const double d1_unslip_upper = d1 - delta;
    const double d1_unslip_lower = d1 - 3 * delta;

    const double d0 = d3 - 3;
    const double d0_upper = d0 + delta;
    const double d0_lower = d0 - delta;
    const double d0_unslip_upper = d0 - delta;
    const double d0_unslip_lower = d0 - 3 * delta;

    for (int i = 0; i < iteration_time; i++) {
      unsliped_config.MoveRelativeDistance(plane_distance_vector);
      // We select four planes, and we move two plan by the burgers vector and
      // calculate the bonds change between the these two planes and the
      // other two.
      std::vector<Atom::Rank> unslipped_atoms_group;
      std::vector<Atom::Rank> slipped_atoms_group;

      for (Atom::Rank j = 0; j < unsliped_config.GetNumAtoms(); ++j) {
        double d_checked =
            DotProduct(unsliped_config.GetAtom(j).relative_position_,
                       plane_index);
        auto check_if_in_between = [](double value, double v_1, double v_2) {
          return (value > std::min(v_1, v_2) && value < std::max(v_1, v_2));
        };
        if (check_if_in_between(d_checked, d0_unslip_lower, d0_unslip_upper) ||
            check_if_in_between(d_checked, d1_unslip_lower, d1_unslip_upper) ||
            check_if_in_between(d_checked, d2_unslip_lower, d2_unslip_upper) ||
            check_if_in_between(d_checked, d3_unslip_lower, d3_unslip_upper)) {
          unslipped_atoms_group.push_back(j);
        }
        if (check_if_in_between(d_checked, d0_lower, d0_upper) ||
            check_if_in_between(d_checked, d1_lower, d1_upper) ||
            check_if_in_between(d_checked, d2_lower, d2_upper) ||
            check_if_in_between(d_checked, d3_lower, d3_upper))
          slipped_atoms_group.push_back(j);
      }
#ifdef MY_DEBUG
      // std::cout << plane_index << std::endl;
      // std::cout << unslipped_atoms_group.size() << std::endl
      //           << slipped_atoms_group.size() << std::endl << std::endl;
#endif
      std::map<Bond, int>
          bonds_map_before = CountBondsBetweenTwoGroupHelper(unsliped_config,
                                                             unslipped_atoms_group,
                                                             slipped_atoms_group);
      //  move atoms
      for (const auto &burgers_vector:burgers_vector_set_) {
        if (DotProduct(plane_index, burgers_vector) != 0)
          continue;
        auto burger_distance_vector = GetBurgerDistanceVector(burgers_vector);
        Config sliped_config = unsliped_config;
        // #ifdef MY_DEBUG
        //         if (i == 0)
        //           unsliped_config.WritePOSCAR("corner_1");
        // #endif
        for (const auto &index:slipped_atoms_group) {
          sliped_config
              .MoveOneAtomRelativeDistance(index, burger_distance_vector);
        }
        // #ifdef MY_DEBUG
        //         if (i == 0)
        //           sliped_config.WritePOSCAR("corner_2");
        // #endif
        std::map<Bond, int>
            bonds_map_after = CountBondsBetweenTwoGroupHelper(sliped_config,
                                                              unslipped_atoms_group,
                                                              slipped_atoms_group);
#ifdef MY_DEBUG
        auto bonds_map_temp = bonds_map_after;
        for (const auto&[key, count] : bonds_map_before) {
          bonds_map_temp[key] -= count;
        }
        for (const auto&[key, count] : bonds_map_temp) {
          std::cout << "#" << key << " " << count << '\n';
        }
        std::cout << '\n';

#endif
        for (const auto&[key, count] : bonds_map_after) {
          bonds_changed[key] += count;
        }
        for (const auto&[key, count] : bonds_map_before) {
          bonds_changed[key] -= count;
        }
      }
    }
  }

#ifdef MY_RELEASE
  for (const auto&[key, count] : bonds_changed) {
    std::cout << "#" << key << " " << count << '\n';
  }
  std::cout << '\n';
#endif
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

double BondCounter::FindD3Helper(const Vector3<double> &plane_index,
                                 const Vector3<double> &box_low_bound,
                                 const Vector3<double> &box_high_bound) {

  return std::max(plane_index.x * box_low_bound.x,
                  plane_index.x * box_high_bound.x)
      + std::max(plane_index.y * box_low_bound.y,
                 plane_index.y * box_high_bound.y)
      + std::max(plane_index.z * box_low_bound.z,
                 plane_index.z * box_high_bound.z);
}

[[maybe_unused]] void BondCounter::GetAtomListBetweenPlanesHelper(
    std::vector<Atom::Rank> &rank_list,
    const Config &config,
    const Vector3<double> &abc,
    const double &d1,
    const double &d2) const {
  double larger = std::max(d1, d2);
  double smaller = std::min(d1, d2);
  for (Atom::Rank i = 0; i < config.GetNumAtoms(); ++i) {
    double d_checked = DotProduct(config.GetAtom(i).relative_position_, abc);
    if (d_checked > larger || d_checked < smaller)
      continue;
    rank_list.push_back(i);
  }
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