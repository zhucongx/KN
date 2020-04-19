#ifndef KN_SRC_VECTORMATRIX_H_
#define KN_SRC_VECTORMATRIX_H_
#include <ostream>
#include <numeric>
#include "armadillo"
// By default, it is always a 1 by 3 vector
template <class NumberType>
struct Vector3 {
  NumberType x, y, z;
  Vector3() = default;;
  Vector3(NumberType x, NumberType y, NumberType z) : x(x), y(y), z(z) {};
  Vector3 &operator+=(const Vector3 &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  };
  Vector3 &operator-=(const Vector3 &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  };
  friend std::ostream &operator<<(std::ostream &os, const Vector3 &vector_3) {
    os << vector_3.x << " " << vector_3.y << " " << vector_3.z;
    return os;
  };
  bool operator<(const Vector3 &rhs) const {
    if (x < rhs.x)
      return true;
    if (rhs.x < x)
      return false;
    if (y < rhs.y)
      return true;
    if (rhs.y < y)
      return false;
    return z < rhs.z;
  };

  [[nodiscard]] Vector3<int> ConvertToInt() const {
    return {(static_cast<int>(x)),
            (static_cast<int>(y)),
            (static_cast<int>(z))};
  };
  [[nodiscard]] Vector3<double> ConvertToDouble() const {
    return {(static_cast<double>(x)),
            (static_cast<double>(y)),
            (static_cast<double>(z))};
  };
};

template <class NumberType>
struct Matrix33 {
  Vector3<NumberType> row1, row2, row3;
};
template <class NumberType>
inline NumberType Max(const Vector3<NumberType> &vector) {
  return std::max(std::max(vector.x, vector.y), vector.z);
}

template <class NumberType>
inline NumberType Min(const Vector3<NumberType> &vector) {
  return std::min(std::min(vector.x, vector.y), vector.z);
}

template <class NumberType>
inline NumberType Sum(const Vector3<NumberType> &vector) {
  return vector.x + vector.y + vector.z;
}

inline Vector3<double> Floor(const Vector3<double> &vector) {
  return {floor(vector.x), floor(vector.y), floor(vector.z)};
}

inline int GCD(const Vector3<int> &vector) {
  return std::gcd(vector.x, std::gcd(vector.y, vector.z));
}

template <class NumberType>
inline bool operator==(const Vector3<NumberType> &lhs,
                       const Vector3<NumberType> &rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}
template <class NumberType>
inline bool operator!=(const Vector3<NumberType> &lhs,
                       const Vector3<NumberType> &rhs) {
  return !(rhs == lhs);
}
template <class NumberType>
inline Vector3<NumberType> operator+(const Vector3<NumberType> &lhs,
                                     const Vector3<NumberType> &rhs) {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}
template <class NumberType>
inline Vector3<NumberType> operator-(const Vector3<NumberType> &lhs,
                                     const Vector3<NumberType> &rhs) {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}
template <class NumberType>
inline Vector3<NumberType> operator*(const Vector3<NumberType> &vector,
                                     const NumberType &factor) {
  return {vector.x * factor, vector.y * factor, vector.z * factor};
}
template <class NumberType>
inline Vector3<NumberType> operator*(const NumberType &factor,
                                     const Vector3<NumberType> &vector) {
  return operator*(vector, factor);
}
template <class NumberType>
inline Vector3<NumberType> CrossProduct(const Vector3<NumberType> &first,
                                        const Vector3<NumberType> &second) {
  return {first.y * second.z - first.z * second.y,
          first.z * second.x - first.x * second.z,
          first.x * second.y - first.y * second.x};
}
template <class NumberType>
inline NumberType DotProduct(const Vector3<NumberType> &first,
                             const Vector3<NumberType> &second) {
  return first.x * second.x + first.y * second.y + first.z * second.z;
}
template <class NumberType>
inline NumberType InnerProduct(const Vector3<NumberType> &vector) {
  return vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
}
template <class NumberType>
inline Vector3<NumberType> StarProduct(const Vector3<NumberType> &first,
                                       const Vector3<NumberType> &second) {
  return {first.x * second.x, first.y * second.y, first.z * second.z};
}
template <class NumberType>
inline Vector3<NumberType> StarDivide(const Vector3<NumberType> &dividend,
                                      const Vector3<NumberType> &divisor) {
  return {dividend.x / divisor.x, dividend.y / divisor.y,
          dividend.z / divisor.z};
}
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
inline Matrix33<NumberType> InverseMatrix33(const Matrix33<NumberType> &input) {
  arma::mat mat_input = {{input.row1.x, input.row1.y, input.row1.z},
                         {input.row2.x, input.row2.y, input.row2.z},
                         {input.row3.x, input.row3.y, input.row3.z}};
  arma::mat inverse_matrix = arma::inv(mat_input);
  return {{inverse_matrix(0, 0), inverse_matrix(0, 1), inverse_matrix(0, 2)},
          {inverse_matrix(1, 0), inverse_matrix(1, 1), inverse_matrix(1, 2)},
          {inverse_matrix(2, 0), inverse_matrix(2, 1), inverse_matrix(2, 2)}};
  //
  // return {{inverse_matrix(0, 0), inverse_matrix(0, 1), inverse_matrix(0, 2)},
  //         {inverse_matrix(1, 0), inverse_matrix(1, 1), inverse_matrix(1, 2)},
  //         {inverse_matrix(2, 0), inverse_matrix(2, 1), inverse_matrix(2, 2)}};

}
#endif //KN_SRC_VECTORMATRIX_H_
