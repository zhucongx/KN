#ifndef KN_INCLUDE_CONFIG_H_
#define KN_INCLUDE_CONFIG_H_

#include <string>
#include <array>
#include <utility>
#include <vector>
#include <map>
#include "Atom.h"
#include "Bond.h"
namespace box {
class Config {
  // Todo output neighbor information
 public:
  Config();
  bool operator<(const Config &rhs) const;
  void Initialize();
  [[nodiscard]] bool IsCubic() const;
  void ConvertRelativeToCartesian();
  void ConvertCartesianToRelative();
  void Perturb();
  virtual void UpdateNeighbors(double first_r_cutoff, double second_r_cutoff);
  // update both atoms' relative and absolute positions according to periodic
  // boundary condition
  void WrapRelativePosition();
  // void WrapAbsolutePosition();
  // void ShiftAtomToCentral(const Atom::Rank &id);

  void MoveRelativeDistance(const Vector3 &distance_vector);
  void MoveOneAtomRelativeDistance(const Atom::Rank &index, const Vector3 &distance_vector);

  // void MoveAbsoluteDistance(const Vector3 &distance_vector);
  std::map<Bond, int> CountAllBonds(double r_cutoff);
  void ReadConfig(const std::string &file_name);
  void ReadPOSCAR(const std::string &file_name);
  void WriteConfig(const std::string &file_name = "config") const;
  // Write Configuration out as POSCAR file. If the show_vacancy_option is
  // true, output will have "X" for visualization. If false, vacancies will be
  // ignored for VASP calculation.
  void WritePOSCAR(const std::string &file_name = "POSCAR",
                   const bool &show_vacancy_option = false) const;

  void GenerateUnitCell(const Matrix33 &bravais_matrix,
                        const std::vector<std::pair<std::string, Vector3>> &type_position_list);
  void Duplicate(const std::array<int, kDimension> &factors);
  void GenerateFCC(const double &lattice_constant_a,
                   const std::string &element,
                   const std::array<int, kDimension> &factors);
  void GenerateBCC(const double &lattice_constant_a,
                   const std::string &element,
                   const std::array<int, kDimension> &factors);
  void GenerateHCP(const double &lattice_constant_a,
                   const double &lattice_constant_c,
                   const std::string &element,
                   const std::array<int, kDimension> &factors);

  void SetScale(double scale);
  [[nodiscard]] double GetScale() const;

  [[nodiscard]] const Matrix33 &GetBasis() const;
  void SetBasis(const Matrix33 &basis);

  [[nodiscard]] const Atom &GetAtom(const int &index) const;
  void AppendAtom(const Atom& atom);
  [[nodiscard]] int GetNumAtoms() const;
  [[nodiscard]] const std::map<std::string, std::vector<Atom::Rank>> &GetElementListMap() const;


 protected:
  double scale_{};
  // double lowx, lowy, lowz, highx, highy, highz, xy xz yz;
  // std::array<double, 9> cell;
  // length of three edges
  // Vector3 length;
  // This three vectors form a matrix matching the matrix in Config
  // and POSCAR file representing three Bravais lattice vectors that can be
  // used to convert relative position to absolute
  Matrix33 basis_{};

  // Three translational Bravais lattice vector
  // Matrix33 reciprocal_matrix_{},

  double energy_{};
  /// The index of atom in the vector is always same as of the id of the atom ?
  std::vector<Atom> atom_list_;
  // indicate if the Config has found Atoms' neighbor list
  bool neighbor_found_{};

  std::map<std::string, std::vector<Atom::Rank>> element_list_map_;
};

}// namespace box

#endif //KN_INCLUDE_CONFIG_H_
