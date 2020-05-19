#ifndef KN_INCLUDE_MATRIX33_H_
#define KN_INCLUDE_MATRIX33_H_

// #define ARMA_ALLOW_FAKE_GCC
// #define ARMA_DONT_USE_WRAPPER
// Uncomment this line if there is a compilation error
// #include <armadillo>

#include "Vector3.h"

typedef std::array<Vector3, kDimension> Matrix33;
inline Matrix33 &operator*=(Matrix33 &lhs, const double &factor)
{
  lhs[kXDimension] *= factor;
  lhs[kYDimension] *= factor;
  lhs[kZDimension] *= factor;
  return lhs;
}
inline Matrix33 &operator/=(Matrix33 &lhs, const double &divisor)
{
  lhs[kXDimension] /= divisor;
  lhs[kYDimension] /= divisor;
  lhs[kZDimension] /= divisor;
  return lhs;
}

inline Matrix33 operator*(const Matrix33 &matrix, const double &factor)
{
  Matrix33 temp(matrix);
  return (temp *= factor);
}

inline Matrix33 operator*(const double &factor, const Matrix33 &matrix)
{
  return operator*(matrix, factor);
}
inline Matrix33 operator/(const Matrix33 &matrix, const double &divisor)
{
  Matrix33 temp(matrix);
  return (temp /= divisor);
}

inline Vector3 operator*(const Vector3 &lhs, const Matrix33 &rhs)
{
  return {lhs[kXDimension] * rhs[kXDimension][kXDimension]
              + lhs[kYDimension] * rhs[kYDimension][kXDimension]
              + lhs[kZDimension] * rhs[kZDimension][kXDimension],
          lhs[kXDimension] * rhs[kXDimension][kYDimension]
              + lhs[kYDimension] * rhs[kYDimension][kYDimension]
              + lhs[kZDimension] * rhs[kZDimension][kYDimension],
          lhs[kXDimension] * rhs[kXDimension][kZDimension]
              + lhs[kYDimension] * rhs[kYDimension][kZDimension]
              + lhs[kZDimension] * rhs[kZDimension][kZDimension]};
}

inline Matrix33 InverseMatrix33(const Matrix33 &input)
{
  // arma::mat mat_input = {{input[kXDimension][kXDimension], input[kXDimension][kYDimension], input[kXDimension][kZDimension]},
  //                        {input[kYDimension][kXDimension], input[kYDimension][kYDimension], input[kYDimension][kZDimension]},
  //                        {input[kZDimension][kXDimension], input[kZDimension][kYDimension], input[kZDimension][kZDimension]}};
  // arma::mat inverse_matrix = arma::inv(mat_input);
  // return {{inverse_matrix(0, 0), inverse_matrix(0, 1), inverse_matrix(0, 2)},
  //         {inverse_matrix(1, 0), inverse_matrix(1, 1), inverse_matrix(1, 2)},
  //         {inverse_matrix(2, 0), inverse_matrix(2, 1), inverse_matrix(2, 2)}};
  double det = (input[kXDimension][kXDimension] * input[kYDimension][kYDimension]
      * input[kZDimension][kZDimension]
      - input[kXDimension][kXDimension] * input[kYDimension][kZDimension]
          * input[kZDimension][kYDimension]
      - input[kXDimension][kYDimension] * input[kYDimension][kXDimension]
          * input[kZDimension][kZDimension]
      + input[kXDimension][kYDimension] * input[kYDimension][kZDimension]
          * input[kZDimension][kXDimension]
      + input[kXDimension][kZDimension] * input[kYDimension][kXDimension]
          * input[kZDimension][kYDimension]
      - input[kXDimension][kZDimension] * input[kYDimension][kYDimension]
          * input[kZDimension][kXDimension]);
  return {
      {{(input[kYDimension][kYDimension] * input[kZDimension][kZDimension]
          - input[kYDimension][kZDimension] * input[kZDimension][kYDimension]) / det,
        (input[kXDimension][kZDimension] * input[kZDimension][kYDimension]
            - input[kXDimension][kYDimension] * input[kZDimension][kZDimension]) / det,
        (input[kXDimension][kYDimension] * input[kYDimension][kZDimension]
            - input[kXDimension][kZDimension] * input[kYDimension][kYDimension]) / det},
       {(input[kYDimension][kZDimension] * input[kZDimension][kXDimension]
           - input[kYDimension][kXDimension] * input[kZDimension][kZDimension]) / det,
        (input[kXDimension][kXDimension] * input[kZDimension][kZDimension]
            - input[kXDimension][kZDimension] * input[kZDimension][kXDimension]) / det,
        (input[kXDimension][kZDimension] * input[kYDimension][kXDimension]
            - input[kXDimension][kXDimension] * input[kYDimension][kZDimension]) / det},
       {(input[kYDimension][kXDimension] * input[kZDimension][kYDimension]
           - input[kYDimension][kYDimension] * input[kZDimension][kXDimension]) / det,
        (input[kXDimension][kYDimension] * input[kZDimension][kXDimension]
            - input[kXDimension][kXDimension] * input[kZDimension][kYDimension]) / det,
        (input[kXDimension][kXDimension] * input[kYDimension][kYDimension]
            - input[kXDimension][kYDimension] * input[kYDimension][kXDimension]) / det}}};
}
#endif //KN_INCLUDE_MATRIX33_H_




