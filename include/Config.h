#ifndef KN_INCLUDE_CONFIG_H_
#define KN_INCLUDE_CONFIG_H_

#include <string>
#include <iostream>
#include <sstream>
#include <array>
#include <utility>
#include <vector>
#include <map>
#include "Atom.h"
#include "Bond.h"
namespace box
{
class Config
{
  // Todo finsih GenerateUnitCell
 public:
  Config();
  bool operator<(const Config &rhs) const;
  void Initialize();
  [[nodiscard]] bool IsCubic() const;
  void ConvertRelativeToAbsolute();
  void ConvertAbsoluteToRelative();
  void Perturb();
  virtual void UpdateNeighbors(double first_r_cutoff, double second_r_cutoff);
  // update both atoms' relative and absolute positions according to periodic
  // boundary condition
  void WrapRelativePosition();
  // void WrapAbsolutePosition();
  // void ShiftAtomToCentral(const Atom::Rank &id);

  void MoveRelativeDistance(const Vector3<double> &distance_vector);
  void MoveOneAtomRelativeDistance(const Atom::Rank &index, const Vector3<double> &distance_vector);

  // void MoveAbsoluteDistance(const Vector3<double> &distance_vector);
  std::map<Bond, int> CountAllBonds(double r_cutoff);
  bool ReadConfig(const std::string &file_name);
  bool ReadPOSCAR(const std::string &file_name);
  void WriteConfig(const std::string &file_name = "config") const;
  // Write Configuration out as POSCAR file. If the show_vacancy_option is
  // true, output will have "X" for visualization. If false, vacancies will be
  // ignored for VASP calculation.
  void WritePOSCAR(const std::string &file_name = "POSCAR",
                   const bool &show_vacancy_option = false) const;

  void GenerateUnitCell(const Matrix33<double> &bravais_matrix,
                        const std::vector<std::pair<std::string, Vector3<double>>> &type_position_list);
  void Duplicate(const Vector3<int> &factors);
  void GenerateFCC(const double &lattice_constant_a, const std::string &element, const Vector3<int> &factors);
  void GenerateBCC(const double &lattice_constant_a, const std::string &element, const Vector3<int> &factors);
  void GenerateHCP(const double &lattice_constant_a,
                   const double &lattice_constant_c,
                   const std::string &element,
                   const Vector3<int> &factors);
  [[nodiscard]] const Matrix33<double> &GetBravaisMatrix() const;
  [[nodiscard]] const Matrix33<double> &GetInverseBravaisMatrix() const;
  [[nodiscard]] const Atom &GetAtom(const Atom::Rank &index) const;
  [[nodiscard]] int GetNumAtoms() const;
 protected:
  double scale_{};
  // double lowx, lowy, lowz, highx, highy, highz, xy xz yz;
  // std::array<double, 9> cell;
  // length of three edges
  // Vector3<double> length;
  // This three vectors form a matrix matching the matrix in Config
  // and POSCAR file representing three Bravais lattice vectors that can be
  // used to convert relative position to absolute
  Matrix33<double> bravais_matrix_{};
  // The inverse three Bravais lattice vectors that can be used to convert
  // absolute to relative position
  Matrix33<double> inverse_bravais_matrix_{};

  // Three translational Bravais lattice vector
  // Matrix33<double> reciprocal_matrix_{},
  int num_atoms_{};
  double energy_{};
  // The index of atom in the vector is always same as of the id of the atom
  std::vector<Atom> atom_list_;
  // indicate if the Config has found Atoms' neighbor list
  bool neighbor_found_{};
  std::map<std::string, std::vector<Atom::Rank>> element_list_set_;
};

}// namespace box

#endif //KN_INCLUDE_CONFIG_H_
