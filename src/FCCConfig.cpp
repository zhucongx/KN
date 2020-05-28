#include "FCCConfig.h"

namespace box {
// todo need to be fixed
void FCCConfig::UpdateNeighbors(double first_r_cutoff,
                                double second_r_cutoff) {
  Vector3 length = {basis_[kXDimension][kXDimension],
                    basis_[kYDimension][kYDimension],
                    basis_[kZDimension][kZDimension]};
  // if the box is a cubic box, we just need to compare relative distance
  const bool cubic_status = IsCubic();
  if (cubic_status) {
    first_r_cutoff /= length[kXDimension];
    second_r_cutoff /= length[kXDimension];
  }
  double first_r_cutoff_square = first_r_cutoff * first_r_cutoff;
  double second_r_cutoff_square = second_r_cutoff * second_r_cutoff;

  for (auto it1 = atom_list_.begin(); it1 < atom_list_.end(); ++it1) {
    for (auto it2 = atom_list_.begin(); it2 < it1; ++it2) {
      Vector3 distance_vector = GetRelativeDistanceVector(*it1, *it2);
      // if the box is not a cubic box, compare absolute distance
      if (!cubic_status) {
        distance_vector = ElementProduct(distance_vector, length);
      }
      if (distance_vector[kXDimension] > second_r_cutoff_square)
        continue;
      if (distance_vector[kYDimension] > second_r_cutoff_square)
        continue;
      if (distance_vector[kZDimension] > second_r_cutoff_square)
        continue;
      double distance_square = Inner(distance_vector);
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
  neighbor_found_ = true;
}

}// namespace box

