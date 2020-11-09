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
#include "Cluster.hpp"
namespace cfg {
class Atom;
class Config {
  public:
    /// Constructor
    Config();
    Config(const Matrix_t &basis, size_t atom_size = 0);
    Config(const Matrix_t &basis, std::vector<Atom> atom_list);
    /// Getter
    [[nodiscard]] size_t GetNumAtoms() const;
    [[nodiscard]] const Matrix_t &GetBasis() const;
    [[nodiscard]] const std::vector<Atom> &GetAtomList() const;
    [[nodiscard]] std::set<std::string> GetTypeSet() const;
    [[nodiscard]] std::map<std::string, std::vector<size_t>> GetElementListMap() const;
    /// Update atoms positions
    void ConvertRelativeToCartesian();
    void ConvertCartesianToRelative();
    void ScaleWith(double scale);
    // update both atoms' relative and cartesian positions according to periodic boundary condition
    void WrapAtomRelative();
    void WrapAtomCartesian();
    // for better performance, shouldn't call Wrap function
    void MoveRelativeDistance(const Vector_t &distance_vector);
    void MoveOneAtomRelativeDistance(size_t index, const Vector_t &distance_vector);
    // add small perturbation to break perfect fcc symmetry this method is about to increase
    // the chance to find lower ground states for VASP software
    void Perturb(std::mt19937_64 &generator);
    void ClearNeighbors();
    // TODO rewrite this function
    void UpdateNeighbors(double first_r_cutoff = Al_const::kFirstNearestNeighborsCutoff,
                         double second_r_cutoff = Al_const::kSecondNearestNeighborsCutoff,
                         double third_r_cutoff = Al_const::kThirdNearestNeighborsCutoff);
    /// Add new atoms
    void AppendAtomWithoutChangingAtomID(const Atom &atom);
    void AppendAtomWithChangingAtomID(Atom atom);
    /// Modify atoms
    void ChangeAtomTypeAt(size_t id, const std::string &type);
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
    double energy_{};
    // The index of atom in the vector is not always same as of the id of the atom
    std::vector<Atom> atom_list_;
  public:
    /// Friend function
    friend void AtomsJump(Config &config, size_t lhs, size_t rhs);
};

// Swap two atoms in a config, and update their near neighbors list
void AtomsJump(Config &config, size_t lhs, size_t rhs);
std::map<Bond, size_t> CountAllBonds(const Config &config);
// Returns the config's type hashmap with the key type name and a categorical label
std::unordered_map<std::string, size_t> GetTypeCategoryHashmap(const Config &config);
std::set<std::string> GetTypeSet(const Config &config);
// jump_pair in a pair of two indexes, this function return the center of these two atoms
Vector_t GetPairCenter(const Config &config, const std::pair<size_t, size_t> &jump_pair);
Matrix_t GetPairRotationMatrix(const Config &config, const std::pair<size_t, size_t> &jump_pair);
void RotateAtomVector(std::vector<Atom> &atom_list, const Matrix_t &rotation_matrix);
// Returns the config with the original IDs.

size_t GetVacancyIndex(const Config &config);
Config GenerateFCC(double lattice_constant_a, const std::string &element, const Factor_t &factors);
}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_CONFIG_H_
