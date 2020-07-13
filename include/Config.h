#ifndef KN_INCLUDE_CONFIG_H_
#define KN_INCLUDE_CONFIG_H_

#include <string>
#include <array>
#include <vector>
#include <map>
#include <random>
#include "Atom.h"
#include "Bond.h"
namespace kn {
class Config {
    /// Todo output neighbor information
  public:

    Config();
    Config(const Matrix33 &basis, int atom_size);
    bool operator<(const Config &rhs) const;
    void ConvertRelativeToCartesian();
    void ConvertCartesianToRelative();
    void UpdateNeighbors(double first_r_cutoff = Al_const::kFirstNearestNeighborsCutoff,
                         double second_r_cutoff = Al_const::kSecondNearestNeighborsCutoff,
                         double third_r_cutoff = Al_const::kThirdNearestNeighborsCutoff);

    // update both atoms' relative and absolute positions according to periodic
    // boundary condition
    void WrapAtomRelative();
    void WrapAtomCartesian();
    void MoveRelativeDistance(const Vector3 &distance_vector);
    void MoveOneAtomRelativeDistance(int index, const Vector3 &distance_vector);
    // add small perturbation to break perfect fcc symmetry this method is about to increase
    // the chance to find lower ground states for VASP software
    void Perturb(std::mt19937_64 &generator);
    [[nodiscard]] int GetNumAtoms() const;

    [[nodiscard]] const Matrix33 &GetBasis() const;

    void AppendAtom(const Atom &atom);
    [[nodiscard]] const std::vector<Atom> &GetAtomList() const;

    [[nodiscard]] const std::map<std::string, std::vector<int>> &GetElementListMap() const;
    [[nodiscard]] bool IsNeighborFound() const;
    void SetNeighborFound(bool neighbor_found);

  protected:
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
    // The index of atom in the vector is not always same as of the id of the atom
    std::vector<Atom> atom_list_;
    // indicate if the Config has found Atoms' neighbor list
    bool neighbor_found_{};

    // using map data structure because we want to keep the order
    std::map<std::string, std::vector<int>> element_list_map_;
};

std::map<Bond, int> CountAllBonds(Config &config);

} // namespace kn
#endif //KN_INCLUDE_CONFIG_H_
