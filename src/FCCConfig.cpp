#include "FCCConfig.h"

namespace box {

void FCCConfig::GenerateFCC(const double &lattice_constant_a,
                            const std::string &element,
                            const Int3 &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  cell_.SetFirstBravaisVector({lattice_constant_a * factors.x, 0, 0});
  cell_.SetSecondBravaisVector({0, lattice_constant_a * factors.y, 0});
  cell_.SetThirdBravaisVector({0, 0, lattice_constant_a * factors.z});
  cell_.SetScale(1.0);
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors.x);
  auto y_length = static_cast<double>(factors.y);
  auto z_length = static_cast<double>(factors.z);
  for (int k = 0; k < factors.z; ++k) {
    for (int j = 0; j < factors.y; ++j) {
      for (int i = 0; i < factors.x; ++i) {
        auto x_reference = static_cast<double>(i);
        auto y_reference = static_cast<double>(j);
        auto z_reference = static_cast<double>(k);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                x_reference / x_length,
                                y_reference / y_length,
                                z_reference / z_length);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                (x_reference + 0.5) / x_length,
                                (y_reference + 0.5) / y_length,
                                z_reference / z_length);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                (x_reference + 0.5) / x_length,
                                y_reference / y_length,
                                (z_reference + 0.5) / z_length);
        atom_list_.emplace_back(num_atoms_++, mass, element,
                                x_reference / x_length,
                                (y_reference + 0.5) / y_length,
                                (z_reference + 0.5) / z_length);
      }
    }
  }
  ConvertRelativeToAbsolute();
}

void FCCConfig::UpdateNeighbors(double first_r_cutoff,
                                double second_r_cutoff) {
  Double3 length = {cell_.GetFirstBravaisVector().x,
                    cell_.GetSecondBravaisVector().y,
                    cell_.GetThirdBravaisVector().z};
  // if the box is a cubic box, we just need to compare relative distance
  const bool cubic_status = cell_.IsCubic();
  if (cubic_status) {
    first_r_cutoff /= length.x;
    second_r_cutoff /= length.x;
  }
  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;

  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Double3 distance_vector = GetRelativeDistanceVector(*it1, *it2);
      // if the box is not a cubic box, compare absolute distance
      if (!cubic_status) {
        distance_vector = double3_calc::StarProduct(distance_vector, length);
      }
      if (distance_vector.x > second_r_cutoff_square)
        continue;
      if (distance_vector.y > second_r_cutoff_square)
        continue;
      if (distance_vector.z > second_r_cutoff_square)
        continue;
      double distance_square =
          double3_calc::InnerProduct(distance_vector);
      if (distance_square < second_r_cutoff_square) {
        if (distance_square < first_r_cutoff_square) {
          it1->first_nearest_neighbor_list_.emplace_back(it2->GetId());
          it2->first_nearest_neighbor_list_.emplace_back(it1->GetId());
        }
        it1->second_nearest_neighbor_list_.emplace_back(it2->GetId());
        it2->second_nearest_neighbor_list_.emplace_back(it1->GetId());
      }
    }
  }
}

}// namespace box

