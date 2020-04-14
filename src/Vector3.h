#ifndef KN_SRC_VECTOR3_H_
#define KN_SRC_VECTOR3_H_

template <class NumberType>
class Vector3 {
 public:
  Vector3() = default;;
  Vector3(NumberType x, NumberType y, NumberType z) : x(x), y(y), z(z) {};
  Vector3 &operator+=(const Vector3 &rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }
  Vector3 &operator-=(const Vector3 &rhs) {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }
  friend bool operator==(const Vector3 &lhs, const Vector3 &rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
  }
  friend bool operator!=(const Vector3 &lhs, const Vector3 &rhs) {
    return !(rhs == lhs);
  }
  friend Vector3 operator+(const Vector3 &lhs, const Vector3 &rhs ){
    return {lhs.x + rhs.x,lhs.y + rhs.y,lhs.z + rhs.z};
  }
  friend Vector3 operator-(const Vector3 &lhs, const Vector3 &rhs ){
    return {lhs.x - rhs.x,lhs.y - rhs.y,lhs.z - rhs.z};
  }
  friend Vector3 operator*(const Vector3 &vector, const NumberType &factor) {
    return {vector.x * factor, vector.y * factor, vector.z * factor};
  }
  friend Vector3 operator*( const NumberType &factor,const Vector3 &vector) {
    return operator*(vector,factor);
  }
  friend Vector3 CrossProduct(const Vector3 &first, const Vector3 &second) {
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
  }
  friend NumberType DotProduct(const Vector3 &first, const Vector3 &second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
  }
  friend NumberType InnerProduct(const Vector3 &vector) {
    return vector.x * vector.x + vector.y * vector.y + vector.z * vector.z;
  }
  friend Vector3 StarProduct(const Vector3 &first, const Vector3 &second) {
    return {first.x * second.x, first.y * second.y, first.z * second.z};
  }
  // first(1*3) * right(3*3)
  friend Vector3 LinearTransform(const Vector3 &first,
                                 const Vector3 &second1,
                                 const Vector3 &second2,
                                 const Vector3 &second3) {
    return {first.x * second1.x + first.y * second2.x + first.z * second3.x,
            first.x * second1.y + first.y * second2.y + first.z * second3.y,
            first.x * second1.z + first.y * second2.z + first.z * second3.z};
  }
  NumberType x, y, z;
};
#endif //KN_SRC_VECTOR3_H_
