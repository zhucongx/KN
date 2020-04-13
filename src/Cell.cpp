#include "Cell.h"

namespace box {

Cell::Cell() = default;

double Cell::GetScale() const {
  return scale_;
}
void Cell::SetScale(double scale) {
  scale_ = scale;
}

const Double3 &Cell::GetFirstBravaisVector() const {
  return first_bravais_vector_;
}
void Cell::SetFirstBravaisVector(const Double3 &first_bravais_vector) {
  first_bravais_vector_ = first_bravais_vector;
}
const Double3 &Cell::GetSecondBravaisVector() const {
  return second_bravais_vector_;
}
void Cell::SetSecondBravaisVector(const Double3 &second_bravais_vector) {
  second_bravais_vector_ = second_bravais_vector;
}
const Double3 &Cell::GetThirdBravaisVector() const {
  return third_bravais_vector_;
}
void Cell::SetThirdBravaisVector(const Double3 &third_bravais_vector) {
  third_bravais_vector_ = third_bravais_vector;
}

void Cell::Initialize() {
  scale_ = 1.0;
  first_bravais_vector_ = {0, 0, 0};
  second_bravais_vector_ = {0, 0, 0};
  third_bravais_vector_ = {0, 0, 0};
}
bool Cell::IsCubic() const {
  return first_bravais_vector_.x == second_bravais_vector_.y &&
      second_bravais_vector_.y == third_bravais_vector_.z &&
      first_bravais_vector_.y == 0 &&
      first_bravais_vector_.z == 0 &&
      second_bravais_vector_.x == 0 &&
      second_bravais_vector_.z == 0 &&
      third_bravais_vector_.x == 0 &&
      third_bravais_vector_.y == 0;
}

}// namespace box
