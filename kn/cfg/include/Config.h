#ifndef KN_KN_CFG_INCLUDE_CONFIG_H_
#define KN_KN_CFG_INCLUDE_CONFIG_H_
#include <string>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include "Atom.h"
#include "Bond.h"
#include "Clusters.hpp"
namespace cfg {
class Atom;
class Config {
  public:
    /// Constructor
    Config();
    Config(const Matrix_t &basis, int atom_size = 0);
    /// Getter
    [[nodiscard]] int GetNumAtoms() const;
    [[nodiscard]] const Matrix_t &GetBasis() const;
    [[nodiscard]] const std::vector<Atom> &GetAtomList() const;
    [[nodiscard]] std::set<std::string> GetTypeSet() const;
    [[nodiscard]] std::map<std::string, std::vector<int>> GetElementListMap() const;
    /// Update atoms positions
    void ConvertRelativeToCartesian();
    void ConvertCartesianToRelative();
    void ScaleWith(double scale);
    // update both atoms' relative and cartesian positions according to periodic boundary condition
    void WrapAtomRelative();
    void WrapAtomCartesian();
    // for better performance, shouldn't call Wrap function
    void MoveRelativeDistance(const Vector_t &distance_vector);
    void MoveOneAtomRelativeDistance(int index, const Vector_t &distance_vector);
    // add small perturbation to break perfect fcc symmetry this method is about to increase
    // the chance to find lower ground states for VASP software
    void Perturb(std::mt19937_64 &generator);
    // TODO rewrite this function
    void UpdateNeighbors(double first_r_cutoff = Al_const::kFirstNearestNeighborsCutoff,
                         double second_r_cutoff = Al_const::kSecondNearestNeighborsCutoff,
                         double third_r_cutoff = Al_const::kThirdNearestNeighborsCutoff);
    /// Add new atoms
    void AppendAtomWithoutChangingAtomID(const Atom &atom);
    void AppendAtomWithChangingAtomID(Atom atom);
    /// Modify atoms
    void ChangeAtomTypeAt(int id, const std::string &type);
    /// IO Todo: rewrite as friend function
    static Config ReadPOSCAR(const std::string &filename, bool update_neighbors = true);
    static Config ReadConfig(const std::string &filename, bool update_neighbors = true);
    // Write Configuration out as POSCAR file. If the show_vacancy_option is
    // true, output will have "X" for visualization. If false, vacancies will be
    // ignored for VASP calculation.
    static void WritePOSCAR(const Config &config,
                            const std::string &filename,
                            bool show_vacancy_option = false);
    static void WriteConfig(const Config &config,
                            const std::string &filename,
                            bool neighbors_info = true);

  private:
    // double lowx, lowy, lowz, highx, highy, highz, xy xz yz;
    // std::array<double, 9> cell;
    // length of three edges
    // Vector_t length;
    // This three vectors form a matrix matching the matrix in Config
    // and POSCAR file representing three Bravais lattice vectors that can be
    // used to convert relative position to absolute
    Matrix_t basis_{};
    // Three translational Bravais lattice vector
    // Matrix_t reciprocal_matrix_{},
    [[maybe_unused]] double energy_{};
    // The index of atom in the vector is not always same as of the id of the atom
    std::vector<Atom> atom_list_;
  public:
    /// Friend function
    friend void AtomsJump(Config &config, int lhs, int rhs);
};

// Swap two atoms in a config, and update their near neighbors list
void AtomsJump(Config &config, int lhs, int rhs);
std::map<Bond, int> CountAllBonds(const Config &config);
std::unordered_map<std::string, int> GetTypeCategoryHashmap(const Config &config);
std::set<std::string> GetTypeSet(const Config &config);
// jump_pair in a pair of two indexes, this function return the center of these two atoms
Vector_t GetPairCenter(const Config &config, const std::pair<int, int> &jump_pair);
Matrix_t GetPairRotationMatrix(const Config &config, const std::pair<int, int> &jump_pair);
void RotateAtomVector(std::vector<Atom> &atom_list, const Matrix_t &rotation_matrix);

Config GenerateFCC(double lattice_constant_a, const std::string &element, const Factor_t &factors);
}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_CONFIG_H_
