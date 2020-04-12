#ifndef KN_SRC_BOX_H_
#define KN_SRC_BOX_H_

#include <array>

#include "constants.h"

class Box {
 public:
  Box();
  Box(const std::array<double, kDimension> &first_bravais_vector,
      const std::array<double, kDimension> &second_bravais_vector,
      const std::array<double, kDimension> &third_bravais_vector);
  [[nodiscard]] const std::array<double,
                                 kDimension> &GetFirstBravaisVector() const;
  double GetScale() const;
  void SetScale(double scale);
  void SetFirstBravaisVector(const std::array<double,
                                              kDimension> &first_bravais_vector);
  [[nodiscard]] const std::array<double,
                                 kDimension> &GetSecondBravaisVector() const;
  void SetSecondBravaisVector(const std::array<double,
                                               kDimension> &second_bravais_vector);
  [[nodiscard]] const std::array<double,
                                 kDimension> &GetThirdBravaisVector() const;
  void SetThirdBravaisVector(const std::array<double,
                                              kDimension> &third_bravais_vector);
  void Initialize();
 private:
  double scale_{};
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  // std::array<double, 9> cell;
  // length of three edges
  // std::array<double, 3> length;
  // This three vectors form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  std::array<double, kDimension> first_bravais_vector_{},
      second_bravais_vector_{}, third_bravais_vector_{};
  // Three translational Bravais lattice vector
  // std::array<double, kDimension> first_reciprocal_vector_{},
  //     second_reciprocal_vector_{}, third_reciprocal_vector_{};

};

#endif //KN_SRC_BOX_H_
