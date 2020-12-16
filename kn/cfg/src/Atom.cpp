#include "Atom.h"

namespace cfg {
Atom::Atom() = default;
// Set both relative and absolute position, but will be corrected later
Atom::Atom(size_t id, double mass, std::string type, double x, double y, double z) :
    id_(id),
    mass_(mass),
    type_(std::move(type)),
    cartesian_position_{x, y, z},
    relative_position_{x, y, z} {
  first_nearest_neighbors_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  second_nearest_neighbors_list_.reserve(Al_const::kNumSecondNearestNeighbors);
  third_nearest_neighbors_list_.reserve(Al_const::kNumThirdNearestNeighbors);
}
Atom::Atom(size_t id, double mass, std::string type, Vector_t position) :
    id_(id),
    mass_(mass),
    type_(std::move(type)),
    cartesian_position_(position),
    relative_position_(position) {
  first_nearest_neighbors_list_.reserve(Al_const::kNumFirstNearestNeighbors);
  second_nearest_neighbors_list_.reserve(Al_const::kNumSecondNearestNeighbors);
  third_nearest_neighbors_list_.reserve(Al_const::kNumThirdNearestNeighbors);
}
const Vector_t &Atom::GetCartesianPosition() const {
  return cartesian_position_;
}
const Vector_t &Atom::GetRelativePosition() const {
  return relative_position_;
}
size_t Atom::GetId() const {
  return id_;
}
double Atom::GetMass() const {
  return mass_;
}
const std::string &Atom::GetType() const {
  return type_;
}
const std::vector<size_t> &Atom::GetFirstNearestNeighborsList() const {
  return first_nearest_neighbors_list_;
}
const std::vector<size_t> &Atom::GetSecondNearestNeighborsList() const {
  return second_nearest_neighbors_list_;
}
const std::vector<size_t> &Atom::GetThirdNearestNeighborsList() const {
  return third_nearest_neighbors_list_;
}
// bool Atom::operator<(const Atom &rhs) const {
//   return id_ < rhs.id_;
// }
void Atom::SetCartesianPosition(const Vector_t &cartesian_position) {
  cartesian_position_ = cartesian_position;
}
void Atom::SetRelativePosition(const Vector_t &relative_position) {
  relative_position_ = relative_position;
}
void Atom::SetId(size_t id) {
  id_ = id;
}
void Atom::SetType(const std::string &type) {
  mass_ = FindMass(type_ = type);
}
void Atom::AppendFirstNearestNeighborsList(size_t index) {
  first_nearest_neighbors_list_.emplace_back(index);
}
void Atom::AppendSecondNearestNeighborsList(size_t index) {
  second_nearest_neighbors_list_.emplace_back(index);
}
void Atom::AppendThirdNearestNeighborsList(size_t index) {
  third_nearest_neighbors_list_.emplace_back(index);
}
void Atom::CleanNeighborsLists() {
  first_nearest_neighbors_list_.clear();
  second_nearest_neighbors_list_.clear();
  third_nearest_neighbors_list_.clear();
}

Vector_t GetRelativeDistanceVector(const Atom &first, const Atom &second) {
  Vector_t relative_distance_vector = second.GetRelativePosition() - first.GetRelativePosition();
  // periodic boundary conditions
  for (const auto kDim : All_Dimensions) {
    while (relative_distance_vector[kDim] >= 0.5)
      relative_distance_vector[kDim] -= 1;
    while (relative_distance_vector[kDim] < -0.5)
      relative_distance_vector[kDim] += 1;
  }
  return relative_distance_vector;
}
double FindMass(const std::string &elem) {
  /* --------------------------------------------------------------------------
  * Atomic weight data taken from:
  * Pure Appl. Chem., Vol. 83, No. 2, pp. 359–396, 2011.
  * Atomic weights of the elements 2009 (IUPAC Technical Report)
  * -------------------------------------------------------------------------- */
  static const std::vector<double> mass_list = {
      0.0,                                                             // Vacancy
      1.00797,     4.0026,      6.939,       9.012182,    10.811,      // H  - B
      12.01115,    14.0067,     15.9994,     18.9984032,  20.17976,    // C  - Ne
      22.989769,   24.30506,    26.9815386,  28.086,      30.973762,   // Na - P
      32.064,      35.453,      39.948,      39.0983,     40.078,      // S  - Ca
      44.955912,   47.867,      50.9415,     51.9961,     54.938045,   // Sc - Mn
      55.845,      58.933195,   58.6934,     63.546,      65.38,       // Fe - Zn
      69.723,      72.63,       74.9216,     78.96,       79.904,      // Ga - Br
      83.798,      85.4678,     87.62,       88.90585,    91.224,      // Kr - Zr
      92.90638,    95.96,       98.9062,     101.07,      102.9055,    // Nb - Rh
      106.42,      107.8682,    112.411,     114.818,     118.71,      // Pd - Sn
      121.76,      127.6,       126.90447,   131.293,     132.9054519, // Sb - Cs
      137.327,     138.90547,   140.116,     140.90765,   144.242,     // Ba - Nd
      147.,        150.36,      151.964,     157.25,      158.92535,   // Pm - Tb
      162.5,       164.93032,   167.259,     168.93421,   173.054,     // Dy - Yb
      174.9668,    178.49,      180.94788,   183.84,      186.207,     // Lu - Re
      190.23,      192.217,     195.084,     196.966569,  200.59,      // Os - Hg
      204.383,     207.2,       208.9804,    209.,        210.,        // Tl - At
      222.,        223.,        226.025,     227.028,     232.03806,   // Rn - Th
      231.03588,   238.02891,   237.048,     244.,        243.,        // Pa - Am
      247.,        247.,        251.,        252.,        257.,        // Cm - Fm
      258.,        259.,        260.,        261.11,      262.11,      // Md - Db
      263.12,      262.12,      265.,        266.,        269.,        // Sg - Ds
      272.,        285.                                                // Rg - Cn
  };
  static const std::vector<std::string> element_list = {
      "X",
      "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
      "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
      "Sc", "Ti", "V", "Cr",  "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
      "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
      "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
      "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
      "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
      "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
      "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
      "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
      "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Os"
  };
  // static const std::vector<int> atomic_num_list = {
  //     0,
  //     1,   2,   3,   4,   5,   6,   7,   8,   9,   10,
  //     11,  12,  13,  14,  15,  16,  17,  18,  19,  20,
  //     21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
  //     31,  32,  33,  34,  35,  36,  37,  38,  39,  40,
  //     41,  42,  43,  44,  45,  46,  47,  48,  49,  50,
  //     51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
  //     61,  62,  63,  64,  65,  66,  67,  68,  69,  70,
  //     71,  72,  73,  74,  75,  76,  77,  78,  79,  80,
  //     81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
  //     91,  92,  93,  94,  95,  96,  97,  98,  99,  100,
  //     101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
  //     111, 112, 113, 114, 115, 116, 117, 118
  // };
  /*this "it" is indeed a iterator*/
  auto it = std::find(element_list.begin(), element_list.end(), elem);
  // If this element is not in the list, return 0
  if (it == element_list.end())
    return 0;
  return mass_list[static_cast<size_t>(std::distance(element_list.begin(), it))];
}
}// namespace cfg
