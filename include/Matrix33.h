#ifndef KN_INCLUDE_MATRIX33_H_
#define KN_INCLUDE_MATRIX33_H_

#include "Vector3.h"
template <class NumberType>
struct Matrix33 {
  Vector3<NumberType> row1, row2, row3;
};

template <class NumberType>
inline bool operator==(const Matrix33<NumberType> &lhs,
                       const Matrix33<NumberType> &rhs) {
  return lhs.row1 == rhs.row1 && lhs.row2 == rhs.row2 && lhs.row3 == rhs.row3;
}
template <class NumberType>
inline bool operator!=(const Matrix33<NumberType> &lhs,
                       const Matrix33<NumberType> &rhs) {
  return !(rhs == lhs);
}
template <class NumberType>
inline Vector3<NumberType> operator*(const Vector3<NumberType> &lhs,
                                     const Matrix33<NumberType> &rhs) {
  return {lhs.x * rhs.row1.x + lhs.y * rhs.row2.x + lhs.z * rhs.row3.x,
          lhs.x * rhs.row1.y + lhs.y * rhs.row2.y + lhs.z * rhs.row3.y,
          lhs.x * rhs.row1.z + lhs.y * rhs.row2.z + lhs.z * rhs.row3.z};
}
template <class NumberType>
inline Matrix33<NumberType> operator*(const Matrix33<NumberType> &matrix,
                                      const NumberType &factor) {
  return {matrix.row1 * factor, matrix.row2 * factor, matrix.row3 * factor};
}
template <class NumberType>
inline Matrix33<NumberType> operator*(const NumberType &factor,
                                      const Matrix33<NumberType> &matrix) {
  return operator*(matrix, factor);
}
inline Matrix33<double> operator/(const Matrix33<double> &matrix,
                                      const double &factor) {
  return {matrix.row1 / factor, matrix.row2 / factor, matrix.row3 / factor};
}
Matrix33<double> InverseMatrix33(const Matrix33<double> &input);

#endif //KN_INCLUDE_MATRIX33_H_
