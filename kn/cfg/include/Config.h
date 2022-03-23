#ifndef KN_KN_CFG_INCLUDE_CONFIG_H_
#define KN_KN_CFG_INCLUDE_CONFIG_H_
#include <string>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include "Atom.h"
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
    [[nodiscard]] const std::unordered_map<size_t, size_t> &GetSiteIdToAtomIdHashmap() const;
    [[nodiscard]] std::set<std::string> GetTypeSet() const;
    [[nodiscard]] std::map<std::string, std::vector<size_t> > GetElementListMap() const;
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
                         double second_r_cutoff = Al_const::kSecondNearestNeighborsCutoffL,
                         double third_r_cutoff = Al_const::kThirdNearestNeighborsCutoff,
                         double fourth_r_cutoff = Al_const::kFourthNearestNeighborsCutoff,
                         double fifth_r_cutoff = Al_const::kFifthNearestNeighborsCutoff,
                         double sixth_r_cutoff = Al_const::kSixthNearestNeighborsCutoff,
                         double seventh_r_cutoff = Al_const::kSeventhNearestNeighborsCutoff);
    void UpdateFirstSecondThirdNeighbors(
        double first_r_cutoff = Al_const::kFirstNearestNeighborsCutoff,
        double second_r_cutoff_l = Al_const::kSecondNearestNeighborsCutoffL,
        double second_r_cutoff_u = Al_const::kSecondNearestNeighborsCutoffU,
        double third_r_cutoff = Al_const::kThirdNearestNeighborsCutoff);
    void CreateSiteIdToAtomIdHashmap();
    /// Add new atoms
    void AppendAtomWithoutChangingAtomID(const Atom &atom);
    void AppendAtomWithChangingAtomID(Atom atom);
    // Removes the atom that has the id
    void RemoveAtomWithID(size_t id);
    /// Modify atoms
    void SortAtomListWith(const std::function<bool(const cfg::Atom &, const cfg::Atom &)> &compare);
    void ReIndexAtoms();
    void ChangeAtomTypeAt(size_t id, const std::string &type);
    /// IO Todo: rewrite as friend function
    static Config ReadPOSCAR(const std::string &filename, bool update_neighbors = true);
    static Config ReadConfig(const std::string &filename, size_t update_neighbors);
    // Write Configuration out as POSCAR file. If the show_vacancy_option is
    // true, output will have "X" for visualization. If false, vacancies will be
    // ignored for VASP calculation.
    static void WritePOSCAR(const Config &config,
                            const std::string &filename,
                            bool show_vacancy_option = false);
    static void WriteConfig(const Config &config,
                            const std::string &filename,
                            size_t neighbors_info);

    // Three translational Bravais lattice vector
    // Matrix_t reciprocal_matrix_{},
    // The index of atom in the vector is not always same as of the id of the atom
    std::vector<Atom> atom_list_{};
  private:
    // double lowx, lowy, lowz, highx, highy, highz, xy xz yz;
    // std::array<double, 9> cell;
    // length of three edges
    // Vector_t length;
    // This three vectors form a matrix matching the matrix in Config
    // and POSCAR file representing three Bravais lattice vectors that can be
    // used to convert relative position to absolute
    Matrix_t basis_{};
    std::unordered_map<size_t, size_t> site_id_to_atom_id_hashmap_{};
  public:
    /// Friend function
    friend void AtomsJump(Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair);
    friend void SitesJump(Config &config, const std::pair<size_t, size_t> &site_id_jump_pair);
};

// Swap two atoms in a config, and update their near neighbors list
std::unordered_set<size_t> GetFirstAndSecondThirdNeighborsSetOfJumpPair(
    const Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair);
// The first indicate to atom
std::map<size_t, size_t> GetAtomIDToSiteIDMapOfFirstThreeNeighborsOfJumpPair(
    const Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair);
size_t GetHashOfAState(const Config &config, size_t vacancy_index);
void AtomsJump(Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair);
void SitesJump(Config &config, const std::pair<size_t, size_t> &site_id_jump_pair);
std::map<std::string, size_t> CountAllType(const Config &config);
// std::map<Bond, size_t> CountAllBonds(const Config &config);
// Returns the config's type hashmap with the key type name and a categorical label
std::unordered_map<std::string, size_t> GetTypeCategoryHashmap(const Config &config);
// std::set<std::string> GetTypeSet(const Config &config);
// atom_id_jump_pair in a pair of two indexes, this function return the center of these two atoms
Vector_t GetPairCenter(const Config &config, const std::pair<size_t, size_t> &atom_id_jump_pair);
Matrix_t GetPairRotationMatrix(const Config &config,
                               const std::pair<size_t, size_t> &atom_id_jump_pair);
void RotateAtomVector(std::vector<Atom> &atom_list, const Matrix_t &rotation_matrix);

size_t GetVacancyIndex(const Config &config);
Config GenerateFCC(double lattice_constant_a, const std::string &element, const Factor_t &factors);
}// namespace cfg
#endif //KN_KN_CFG_INCLUDE_CONFIG_H_
