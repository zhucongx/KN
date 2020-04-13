#ifndef KN_SRC_ATOM_H_
#define KN_SRC_ATOM_H_

#include <fstream>
#include <array>
#include <vector>

#include "constants.h"

class Atom {
 public:
  Atom(int id);
  Atom(int id, double mass, std::string type);
  // Set both relative and absolute position, but will be corrected later
  Atom(int id, double mass, std::string type, double x, double y, double z);
  virtual ~Atom();
  void SetId(int id);
  [[nodiscard]] int GetId() const;
  void SetMass(double mass);
  [[nodiscard]] double GetMass() const;
  void SetType(const std::string &type);
  [[nodiscard]] const std::string &GetType() const;
  void SetAbsolutePosition(const Double3 &absolute_position);
  [[nodiscard]] const Double3 &GetAbsolutePosition() const;
  void SetRelativePosition(const Double3 &relative_position);
  [[nodiscard]] const Double3 &GetRelativePosition() const;

  // First and second nearest neighbor list
  std::vector<Atom> second_near_neighbor_list_;
  // First nearest neighbor list
  std::vector<Atom> first_nearest_neighbor_list_;
 private:
  // atom id which is an unique int for every atom indexed form 0
  int id_;
  double mass_;
  std::string type_;
  // absolute position
  Double3 absolute_position_{};
  // relative position in the box
  Double3 relative_position_{};
};
#endif //KN_SRC_ATOM_H_
