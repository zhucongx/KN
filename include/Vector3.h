#ifndef KN_INCLUDE_VECTOR3_H_
#define KN_INCLUDE_VECTOR3_H_

#include <cmath>

#include <ostream>
#include <numeric>

#include "Constants.h"

// Todo: rewrite using Eigen
// By default, it is always a 1 by 3 vector

typedef std::array<double, kDimension> Vector3;

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

inline static double Max(const Vector3 &vector)
{
  return std::max(std::max(vector[kXDimension], vector[kYDimension]), vector[kZDimension]);
}
inline static double Min(const Vector3 &vector)
{
  return std::min(std::min(vector[kXDimension], vector[kYDimension]), vector[kZDimension]);
}
inline static double Sum(const Vector3 &vector)
{
  return vector[kXDimension] + vector[kYDimension] + vector[kZDimension];
}
inline static Vector3 Abs(const Vector3 &vector)
{
  return {std::abs(vector[kXDimension]), std::abs(vector[kYDimension]),
          std::abs(vector[kZDimension])};
}
inline static Vector3 Floor(const Vector3 &vector)
{
  return {floor(vector[kXDimension]), floor(vector[kYDimension]), floor(vector[kZDimension])};
}
inline static Vector3 CrossProduct(const Vector3 &first, const Vector3 &second)
{
  return {first[kYDimension] * second[kZDimension] - first[kZDimension] * second[kYDimension],
          first[kZDimension] * second[kXDimension] - first[kXDimension] * second[kZDimension],
          first[kXDimension] * second[kYDimension] - first[kYDimension] * second[kXDimension]};
}
inline static double DotProduct(const Vector3 &first, const Vector3 &second)
{
  return first[kXDimension] * second[kXDimension] + first[kYDimension] * second[kYDimension]
      + first[kZDimension] * second[kZDimension];
}
inline static double InnerProduct(const Vector3 &vector)
{
  return vector[kXDimension] * vector[kXDimension] + vector[kYDimension] * vector[kYDimension]
      + vector[kZDimension] * vector[kZDimension];
}
inline static Vector3 StarProduct(const Vector3 &first, const Vector3 &second)
{
  return {first[kXDimension] * second[kXDimension], first[kYDimension] * second[kYDimension],
          first[kZDimension] * second[kZDimension]};
}
inline static Vector3 StarDivide(const Vector3 &dividend, const Vector3 &divisor)
{
  return {dividend[kXDimension] / divisor[kXDimension],
          dividend[kYDimension] / divisor[kYDimension],
          dividend[kZDimension] / divisor[kZDimension]};
}

#endif //KN_INCLUDE_VECTOR3_H_
