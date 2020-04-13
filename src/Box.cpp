#include "Box.h"

Box::Box() = default;

double Box::GetScale() const {
  return scale_;
}
void Box::SetScale(double scale) {
  scale_ = scale;
}

const Double3 &Box::GetFirstBravaisVector() const {
  return first_bravais_vector_;
}
void Box::SetFirstBravaisVector(const Double3 &first_bravais_vector) {
  first_bravais_vector_ = first_bravais_vector;
}
const Double3 &Box::GetSecondBravaisVector() const {
  return second_bravais_vector_;
}
void Box::SetSecondBravaisVector(const Double3 &second_bravais_vector) {
  second_bravais_vector_ = second_bravais_vector;
}
const Double3 &Box::GetThirdBravaisVector() const {
  return third_bravais_vector_;
}
void Box::SetThirdBravaisVector(const Double3 &third_bravais_vector) {
  third_bravais_vector_ = third_bravais_vector;
}

void Box::Initialize() {
  scale_ = 1.0;
  first_bravais_vector_.fill(0);
  second_bravais_vector_.fill(0);
  third_bravais_vector_.fill(0);
}
bool Box::IsCubic() const {
  return first_bravais_vector_[kXDim] == second_bravais_vector_[kYDim] &&
      second_bravais_vector_[kYDim] == third_bravais_vector_[kZDim] &&
      first_bravais_vector_[kYDim] == 0 &&
      first_bravais_vector_[kZDim] == 0 &&
      second_bravais_vector_[kXDim] == 0 &&
      second_bravais_vector_[kZDim] == 0 &&
      third_bravais_vector_[kXDim] == 0 &&
      third_bravais_vector_[kYDim] == 0;
}
