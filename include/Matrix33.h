#ifndef KN_INCLUDE_MATRIX33_H_
#define KN_INCLUDE_MATRIX33_H_

// #define ARMA_ALLOW_FAKE_GCC
// #define ARMA_DONT_USE_WRAPPER
// Uncomment this line if there is a compilation error
// #include <armadillo>

#include "Vector3.h"
template <class NumberType>
struct Matrix33
{
  Vector3<NumberType> row1, row2, row3;
};

template <class NumberType>
inline bool operator==(const Matrix33<NumberType> &lhs,
                       const Matrix33<NumberType> &rhs)
{
  return lhs.row1 == rhs.row1 && lhs.row2 == rhs.row2 && lhs.row3 == rhs.row3;
}
template <class NumberType>
inline bool operator!=(const Matrix33<NumberType> &lhs,
                       const Matrix33<NumberType> &rhs)
{
  return !(rhs == lhs);
}
template <class NumberType>
inline Vector3<NumberType> operator*(const Vector3<NumberType> &lhs,
                                     const Matrix33<NumberType> &rhs)
{
  return {lhs.x * rhs.row1.x + lhs.y * rhs.row2.x + lhs.z * rhs.row3.x,
          lhs.x * rhs.row1.y + lhs.y * rhs.row2.y + lhs.z * rhs.row3.y,
          lhs.x * rhs.row1.z + lhs.y * rhs.row2.z + lhs.z * rhs.row3.z};
}

template <class NumberType>
inline Matrix33<NumberType> operator*(const Matrix33<NumberType> &lhs,
                                      const Matrix33<NumberType> &rhs)
{
  return {lhs.row1 * rhs, lhs.row2 * rhs, lhs.row3 * rhs};
}

template <class NumberType>
inline Matrix33<NumberType> operator*(const Matrix33<NumberType> &matrix,
                                      const NumberType &factor)
{
  return {matrix.row1 * factor, matrix.row2 * factor, matrix.row3 * factor};
}
template <class NumberType>
inline Matrix33<NumberType> operator*(const NumberType &factor,
                                      const Matrix33<NumberType> &matrix)
{
  return operator*(matrix, factor);
}
inline Matrix33<double> operator/(const Matrix33<double> &matrix,
                                  const double &factor)
{
  return {matrix.row1 / factor, matrix.row2 / factor, matrix.row3 / factor};
}
inline Matrix33<double> InverseMatrix33(const Matrix33<double> &input)
{
  // arma::mat mat_input = {{input.row1.x, input.row1.y, input.row1.z},
  //                        {input.row2.x, input.row2.y, input.row2.z},
  //                        {input.row3.x, input.row3.y, input.row3.z}};
  // arma::mat inverse_matrix = arma::inv(mat_input);
  // return {{inverse_matrix(0, 0), inverse_matrix(0, 1), inverse_matrix(0, 2)},
  //         {inverse_matrix(1, 0), inverse_matrix(1, 1), inverse_matrix(1, 2)},
  //         {inverse_matrix(2, 0), inverse_matrix(2, 1), inverse_matrix(2, 2)}};
  double det = (input.row1.x * input.row2.y * input.row3.z
      - input.row1.x * input.row2.z * input.row3.y
      - input.row1.y * input.row2.x * input.row3.z
      + input.row1.y * input.row2.z * input.row3.x
      + input.row1.z * input.row2.x * input.row3.y
      - input.row1.z * input.row2.y * input.row3.x);
  return {{(input.row2.y * input.row3.z - input.row2.z * input.row3.y) / det,
           (input.row1.z * input.row3.y - input.row1.y * input.row3.z) / det,
           (input.row1.y * input.row2.z - input.row1.z * input.row2.y) / det},
          {(input.row2.z * input.row3.x - input.row2.x * input.row3.z) / det,
           (input.row1.x * input.row3.z - input.row1.z * input.row3.x) / det,
           (input.row1.z * input.row2.x - input.row1.x * input.row2.z) / det},
          {(input.row2.x * input.row3.y - input.row2.y * input.row3.x) / det,
           (input.row1.y * input.row3.x - input.row1.x * input.row3.y) / det,
           (input.row1.x * input.row2.y - input.row1.y * input.row2.x) / det}};

}

#endif //KN_INCLUDE_MATRIX33_H_
