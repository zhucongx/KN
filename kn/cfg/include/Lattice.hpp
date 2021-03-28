#ifndef KN_KN_CFG_INCLUDE_LATTICE_H_
#define KN_KN_CFG_INCLUDE_LATTICE_H_
namespace cfg {
struct Lattice {
  public:
    Lattice() = default;
    Lattice(size_t id, double x, double y, double z) {
      id_ = id;
      cartesian_position_ = {x, y, z};
      relative_position_ = {x, y, z};
    }
    Lattice(size_t id, Vector_t position) {
      id_ = id;
      cartesian_position_ = position;
      relative_position_ = position;
    }
    // lattice id
    size_t id_{};
    // absolute position
    Vector_t cartesian_position_{};
    // relative position in the box
    Vector_t relative_position_{};
};
}// namespace cfg


#endif //KN_KN_CFG_INCLUDE_LATTICE_H_
