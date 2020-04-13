#ifndef KN_SRC_CELL_H_
#define KN_SRC_CELL_H_

#include <array>

#include "constants.h"

namespace box {

class Cell {
 public:
  Cell();
  [[nodiscard]] double GetScale() const;
  void SetScale(double scale);
  [[nodiscard]] const Double3 &GetFirstBravaisVector() const;
  void SetFirstBravaisVector(const Double3 &first_bravais_vector);
  [[nodiscard]] const Double3 &GetSecondBravaisVector() const;
  void SetSecondBravaisVector(const Double3 &second_bravais_vector);
  [[nodiscard]] const Double3 &GetThirdBravaisVector() const;
  void SetThirdBravaisVector(const Double3 &third_bravais_vector);
  void Initialize();
  [[nodiscard]] bool IsCubic() const;
 private:
  double scale_{};
  // lowx, lowy, lowz, highx, highy, highz, xy xz yz
  // std::array<double, 9> cell;
  // length of three edges
  // Double3 length;
  // This three vectors form a matrix matching the matrix in Config and POSCAR file
  // representing three Bravais lattice vector
  Double3 first_bravais_vector_{}, second_bravais_vector_{},
      third_bravais_vector_{};
  // Three translational Bravais lattice vector
  // DoubleVecfirst_reciprocal_vector_{},
  //     second_reciprocal_vector_{}, third_reciprocal_vector_{};

};

}// namespace box

#endif //KN_SRC_CELL_H_
