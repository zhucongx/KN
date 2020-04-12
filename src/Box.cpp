#include "Box.h"

Box::Box() = default;

double Box::GetScale() const {
  return scale_;
}
void Box::SetScale(double scale) {
  scale_ = scale;
}

const std::array<double, kDimension> &Box::GetFirstBravaisVector() const {
  return first_bravais_vector_;
}
void Box::SetFirstBravaisVector(const std::array<double,
                                                 kDimension> &first_bravais_vector) {
  first_bravais_vector_ = first_bravais_vector;
}
const std::array<double, kDimension> &Box::GetSecondBravaisVector() const {
  return second_bravais_vector_;
}
void Box::SetSecondBravaisVector(const std::array<double,
                                                  kDimension> &second_bravais_vector) {
  second_bravais_vector_ = second_bravais_vector;
}
const std::array<double, kDimension> &Box::GetThirdBravaisVector() const {
  return third_bravais_vector_;
}
void Box::SetThirdBravaisVector(const std::array<double,
                                                 kDimension> &third_bravais_vector) {
  third_bravais_vector_ = third_bravais_vector;
}

void Box::Initialize() {
  scale_ = 1.0;
  first_bravais_vector_.fill(0);
  second_bravais_vector_.fill(0);
  third_bravais_vector_.fill(0);
}
