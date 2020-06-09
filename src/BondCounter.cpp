#include "BondCounter.h"
#include "AtomUtility.h"
#include "ConfigIO.h"

namespace kn {
BondCounter::BondCounter(const std::string &filename,
                         Vector3 factor,
                         Vector3 plane,
                         Vector3 burger_vector
    ) : factor_(factor) {
  config_ = ConfigIO::ReadConfig(filename);
  SetPlane(plane);
  SetBurgersVector(burger_vector);
}

void BondCounter::SetFactor(const Vector3 &factor) {
  factor_ = factor;
}

void BondCounter::SetPlane(Vector3 miller_index) {
  for (double i : {-1.0, 1.0}) {
    for (double j : {-1.0, 1.0}) {
      for (double k : {-1.0, 1.0}) {
        for (int l = 0; l < 3; l++) {
          miller_index =
              {miller_index[kYDimension], miller_index[kZDimension], miller_index[kXDimension]};
          plane_set_.insert(ElementProduct(miller_index, {i, j, k}));
        }
      }
    }
  }
}

void BondCounter::SetBurgersVector(Vector3 miller_index) {
  for (int l = 0; l < 3; l++) {
    miller_index =
        {miller_index[kYDimension], miller_index[kZDimension], miller_index[kXDimension]};
    for (double i : {-1.0, 1.0}) {
      for (double j : {-1.0, 1.0}) {
        for (double k : {-1.0, 1.0}) {
          burgers_vector_set_.insert(ElementProduct(miller_index, {i, j, k}));
        }
      }
    }
  }
}

std::map<Bond, int> BondCounter::GetBondChange() const {
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
  // int iii = -1;

  std::set<std::pair<int, Vector3>> slip_direction_set;
  for (const auto &plane_index : plane_set_) {
    auto unslipped_config = config_;
    Vector3 plane_distance_vector = GetPlaneDistanceVector(plane_index);
    /// {1/90, 1/90, 1/90}

    const double delta = 0.5 * Sum(ElementAbs(plane_distance_vector));
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
      unslipped_config.MoveRelativeDistance(plane_distance_vector);
      // unslipped_config.WritePOSCAR((std::to_string(++iii) + ".poscar"));
      // int jjj = 0;
      // We select four planes, and we move two plan by the burgers vector and
      // calculate the bonds change between the these two planes and the
      // other two.
      std::vector<int> unslipped_atoms_group;
      std::vector<int> slipped_atoms_group;

      auto upslipped_atom_reference = unslipped_config.GetAtomList();
      for (int j = 0; j < unslipped_config.GetNumAtoms(); ++j) {
        double d_checked =
            Dot(upslipped_atom_reference[j].relative_position_, plane_index);
        auto check_if_in_between = [](double value, double v_1, double v_2) {
          return (value > std::min(v_1, v_2) && value < std::max(v_1, v_2));
        };
        if (check_if_in_between(d_checked, d0_unslip_lower, d0_unslip_upper) ||
          check_if_in_between(d_checked, d1_unslip_lower, d1_unslip_upper) ||
          check_if_in_between(d_checked, d2_unslip_lower, d2_unslip_upper) ||
          check_if_in_between(d_checked, d3_unslip_lower, d3_unslip_upper)) {
          unslipped_atoms_group.push_back(j);
        } else if (check_if_in_between(d_checked, d0_lower, d0_upper) ||
          check_if_in_between(d_checked, d1_lower, d1_upper) ||
          check_if_in_between(d_checked, d2_lower, d2_upper) ||
          check_if_in_between(d_checked, d3_lower, d3_upper)) {
          slipped_atoms_group.push_back(j);
        }
      }
#ifdef MY_DEBUG
      // std::cout << plane_index << std::endl;
      // std::cout << unslipped_atoms_group.size() << std::endl
      //           << slipped_atoms_group.size() << std::endl << std::endl;
#endif
      std::map<Bond, int>
          bonds_map_before = CountBondsBetweenTwoGroupHelper(unslipped_config,
                                                             unslipped_atoms_group,
                                                             slipped_atoms_group);
      //  move atoms
      for (const auto &burgers_vector : burgers_vector_set_) {
        if (Dot(plane_index, burgers_vector) != 0)
          continue;
        auto slip_direction = Cross(plane_index, burgers_vector);

        std::pair<int, Vector3>
            slip_direction_check = std::make_pair(i, slip_direction);

        auto it = slip_direction_set.find(slip_direction_check);
        if (it == slip_direction_set.end()) {
          slip_direction_set.insert(slip_direction_check);
        } else {
          continue;
        }

        auto burger_distance_vector = GetBurgerDistanceVector(burgers_vector);
        Config sliped_config = unslipped_config;

        for (const auto &index : slipped_atoms_group) {
          sliped_config
              .MoveOneAtomRelativeDistance(index, burger_distance_vector);
        }
        // sliped_config
        //     .WritePOSCAR((std::to_string(iii) + "_" + std::to_string(jjj++)
        //         + ".poscar"));

        std::map<Bond, int>
            bonds_map_after = CountBondsBetweenTwoGroupHelper(sliped_config,
                                                              unslipped_atoms_group,
                                                              slipped_atoms_group);
        if (bonds_map_after.empty()) {
          continue;
        }

        for (const auto &[key, count] : bonds_map_after) {
          bonds_changed[key] += count;
        }
        for (const auto &[key, count] : bonds_map_before) {
          bonds_changed[key] -= count;
        }
      }
    }
  }
  return bonds_changed;
}

// relative distance
Vector3 BondCounter::GetPlaneDistanceVector(const Vector3 &plane_index) const {
  double inner = Inner(plane_index);
  Vector3 distance_index = (1.0 / inner * plane_index);
  return ElementDivide(distance_index, factor_);
}

Vector3 BondCounter::GetBurgerDistanceVector(const Vector3 &burger_vector) const {
  return ElementDivide(burger_vector, factor_);
}

double BondCounter::FindD3Helper(const Vector3 &plane_index,
                                 const Vector3 &box_low_bound,
                                 const Vector3 &box_high_bound) {

  return std::max(plane_index[kXDimension] * box_low_bound[kXDimension],
                  plane_index[kXDimension] * box_high_bound[kXDimension])
      + std::max(plane_index[kYDimension] * box_low_bound[kYDimension],
                 plane_index[kYDimension] * box_high_bound[kYDimension])
      + std::max(plane_index[kZDimension] * box_low_bound[kZDimension],
                 plane_index[kZDimension] * box_high_bound[kZDimension]);
}

std::map<Bond, int> BondCounter::CountBondsBetweenTwoGroupHelper(
    const Config &config,
    const std::vector<int> &group1,
    const std::vector<int> &group2) {
  std::map<Bond, int> map_out;
  auto atoms_list_reference = config.GetAtomList();

  for (const auto &index1 : group1) {
    for (const auto &index2 : group2) {
      Vector3 relative_distance_vector =
          AtomUtility::GetRelativeDistanceVector(atoms_list_reference[index1],
                                                 atoms_list_reference[index2]);
      double absolute_distance_squared =
          (Inner(relative_distance_vector * config.GetBasis()));

      if (absolute_distance_squared > 8 && absolute_distance_squared < 9) {
        // std::cout << absolute_distance_squared<<'\n';
        map_out[{
            atoms_list_reference[index1].type_,
            atoms_list_reference[index2].type_
        }]++;
      }
    }
  }
  return map_out;
}

} // namespace kn
