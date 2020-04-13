#include "FCCConfig.h"
void FCCConfig::GenerateFCC(const double &lattice_constant_a,
                            const std::string &element,
                            const Int3 &factors) {
  Initialize();
  double mass = elem_info::FindMass(element);
  box_.SetFirstBravaisVector({lattice_constant_a * factors[kXDim], 0, 0});
  box_.SetSecondBravaisVector({0, lattice_constant_a * factors[kYDim], 0});
  box_.SetThirdBravaisVector({0, 0, lattice_constant_a * factors[kZDim]});
  box_.SetScale(1.0);
  num_atoms_ = 0;
  auto x_length = static_cast<double>(factors[kXDim]);
  auto y_length = static_cast<double>(factors[kYDim]);
  auto z_length = static_cast<double>(factors[kZDim]);
  for (int k = 0; k < factors[kZDim]; ++k) {
    for (int j = 0; j < factors[kYDim]; ++j) {
      for (int i = 0; i < factors[kXDim]; ++i) {
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
  Double3 length = {box_.GetFirstBravaisVector()[kXDim],
                    box_.GetSecondBravaisVector()[kYDim],
                    box_.GetThirdBravaisVector()[kZDim]};
  // if the box is a cubic box, we just need to compare relative distance
  const bool cubic_status = box_.IsCubic();
  if (cubic_status) {
    first_r_cutoff /= length[kXDim];
    second_r_cutoff /= length[kXDim];
  }
  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;

  for (int i = 0; i < num_atoms_; i++) {
    for (int j = 0; j < i; j++) {
      Double3 distance_vector = GetRelativeDistanceVector(i, j);
      // if the box is not a cubic box, compare absolute distance
      if (!cubic_status) {
        distance_vector =
            double3_calc::StarProduct(GetRelativeDistanceVector(i, j), length);
      }
      if (distance_vector[kXDim] > second_r_cutoff_square)
        continue;
      if (distance_vector[kYDim] > second_r_cutoff_square)
        continue;
      if (distance_vector[kZDim] > second_r_cutoff_square)
        continue;
      double distance_square =
          double3_calc::InnerProduct(distance_vector);
      if (distance_square < second_r_cutoff_square) {
        if (distance_square < first_r_cutoff_square) {
          atom_list_[i].first_nearest_neighbor_list_.emplace_back(j);
          atom_list_[j].first_nearest_neighbor_list_.emplace_back(i);
        }
        atom_list_[i].second_near_neighbor_list_.emplace_back(j);
        atom_list_[j].second_near_neighbor_list_.emplace_back(i);
      }
    }
  }
}

