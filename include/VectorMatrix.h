#ifndef KN_INCLUDE_VECTORMATRIX_H_
#define KN_INCLUDE_VECTORMATRIX_H_

#include <cmath>

#include <iostream>
#include <numeric>
// #define ARMA_ALLOW_FAKE_GCC
// #define ARMA_DONT_USE_WRAPPER
// Uncomment this line if there is a compilation error
// #include <armadillo>

const int kDimension = 3;
enum
{
  kXDimension, kYDimension, kZDimension
};

// By default, it is always a 1 by 3 vector

typedef std::array<double, kDimension> Vector3;
typedef std::array<Vector3, kDimension> Matrix33;

inline std::ostream &operator<<(std::ostream &os, const Vector3 &vector)
{
  os << vector[kXDimension] << ' ' << vector[kYDimension] << ' ' << vector[kZDimension];
  return os;
}

inline std::istream &operator>>(std::istream &is, Vector3 &vector)
{
  is >> vector[kXDimension] >> vector[kYDimension] >> vector[kZDimension];
  return is;
}

inline bool operator==(const Vector3 &lhs, const Vector3 &rhs)
{
  return lhs[kXDimension] == rhs[kXDimension] &&
      lhs[kYDimension] == rhs[kYDimension] &&
      lhs[kZDimension] == rhs[kZDimension];
}

inline bool operator!=(const Vector3 &lhs, const Vector3 &rhs)
{
  return !(rhs == lhs);
}

inline bool operator<(const Vector3 &lhs, const Vector3 &rhs)
{
  if (lhs[kXDimension] < rhs[kXDimension])
    return true;
  if (rhs[kXDimension] < lhs[kXDimension])
    return false;
  if (lhs[kYDimension] < rhs[kYDimension])
    return true;
  if (rhs[kYDimension] < lhs[kYDimension])
    return false;
  return lhs[kZDimension] < rhs[kZDimension];
}

inline Vector3 &operator+=(Vector3 &lhs, const Vector3 &rhs)
{
  lhs[kXDimension] += rhs[kXDimension];
  lhs[kYDimension] += rhs[kYDimension];
  lhs[kZDimension] += rhs[kZDimension];
  return lhs;
}
inline Vector3 &operator-=(Vector3 &lhs, const Vector3 &rhs)
{
  lhs[kXDimension] -= rhs[kXDimension];
  lhs[kYDimension] -= rhs[kYDimension];
  lhs[kZDimension] -= rhs[kZDimension];
  return lhs;
}
inline Vector3 &operator*=(Vector3 &lhs, const double &factor)
{
  lhs[kXDimension] *= factor;
  lhs[kYDimension] *= factor;
  lhs[kZDimension] *= factor;
  return lhs;
}
inline Vector3 &operator/=(Vector3 &lhs, const double &divisor)
{
  lhs[kXDimension] /= divisor;
  lhs[kYDimension] /= divisor;
  lhs[kZDimension] /= divisor;
  return lhs;
}
inline Vector3 operator+(const Vector3 &lhs, const Vector3 &rhs)
{
  Vector3 temp(lhs);
  return (temp += rhs);
}

inline Vector3 operator-(const Vector3 &lhs, const Vector3 &rhs)
{
  Vector3 temp(lhs);
  return (temp -= rhs);
}

inline Vector3 operator*(const Vector3 &vector, const double &factor)
{
  Vector3 temp(vector);
  return (temp *= factor);
}

inline Vector3 operator*(const double &factor, const Vector3 &vector)
{
  return operator*(vector, factor);
}

inline Vector3 operator/(const Vector3 &vector, const double &divisor)
{
  Vector3 temp(vector);
  return (temp /= divisor);
}

inline double Max(const Vector3 &vector)
{
  return std::max(std::max(vector[kXDimension], vector[kYDimension]), vector[kZDimension]);
}
inline double Min(const Vector3 &vector)
{
  return std::min(std::min(vector[kXDimension], vector[kYDimension]), vector[kZDimension]);
}
inline double Sum(const Vector3 &vector)
{
  return vector[kXDimension] + vector[kYDimension] + vector[kZDimension];
}
inline Vector3 ElementAbs(const Vector3 &vector)
{
  return {std::abs(vector[kXDimension]), std::abs(vector[kYDimension]),
          std::abs(vector[kZDimension])};
}
inline Vector3 ElementFloor(const Vector3 &vector)
{
  return {floor(vector[kXDimension]), floor(vector[kYDimension]), floor(vector[kZDimension])};
}
inline Vector3 Cross(const Vector3 &first, const Vector3 &second)
{
  return {first[kYDimension] * second[kZDimension] - first[kZDimension] * second[kYDimension],
          first[kZDimension] * second[kXDimension] - first[kXDimension] * second[kZDimension],
          first[kXDimension] * second[kYDimension] - first[kYDimension] * second[kXDimension]};
}
inline double Dot(const Vector3 &first, const Vector3 &second)
{
  return first[kXDimension] * second[kXDimension] + first[kYDimension] * second[kYDimension]
      + first[kZDimension] * second[kZDimension];
}
inline double Inner(const Vector3 &vector)
{
  return vector[kXDimension] * vector[kXDimension] + vector[kYDimension] * vector[kYDimension]
      + vector[kZDimension] * vector[kZDimension];
}
inline Vector3 ElementProduct(const Vector3 &first, const Vector3 &second)
{
  return {first[kXDimension] * second[kXDimension], first[kYDimension] * second[kYDimension],
          first[kZDimension] * second[kZDimension]};
}
inline Vector3 ElementDivide(const Vector3 &dividend, const Vector3 &divisor)
{
  return {dividend[kXDimension] / divisor[kXDimension],
          dividend[kYDimension] / divisor[kYDimension],
          dividend[kZDimension] / divisor[kZDimension]};
}
inline double ScalarLength(const Vector3 &vector)
{
  return std::sqrt(Inner(vector));
}
inline Vector3 Normalize(const Vector3 &vector)
{
  double factor = 1.0 / ScalarLength(vector);
  return vector * factor;
}

inline std::ostream &operator<<(std::ostream &os, const Matrix33 &matrix)
{
  os << matrix[kXDimension] << '\n' << matrix[kYDimension] << '\n' << matrix[kZDimension];
  return os;
}

inline std::istream &operator>>(std::istream &is, Matrix33 &matrix)
{
  is >> matrix[kXDimension] >> matrix[kYDimension] >> matrix[kZDimension];
  return is;
}

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
#endif //KN_INCLUDE_VECTORMATRIX_H_
