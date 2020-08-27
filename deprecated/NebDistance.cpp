#include "NebDistance.h"

namespace kn::NebDistance {
static std::vector<std::pair<int, double>> GetSortedIndexDistanceSquaredHelper(
    const Config &config1,
    const Config &config2) {
  std::vector<std::pair<int, double>> index_and_distance_squared_vector;
  for (const auto &atom1 : config1.GetAtomList()) {
    int index = atom1.GetId();
    const auto &atom2 = config2.GetAtomList()[index];
    index_and_distance_squared_vector.emplace_back(index,
                                                   Inner(GetRelativeDistanceVector(atom1, atom2)
                                                             * config1.GetBasis()));
  }

  std::sort(index_and_distance_squared_vector.begin(), index_and_distance_squared_vector.end(),
            [](const auto &first_row, const auto &second_row) {
              return first_row.second < second_row.second;
            });

  return index_and_distance_squared_vector;
}

void PrintTheDistanceFromTwoPOSCARFiles(const std::string &filename1,
                                        const std::string &filename2) {
  auto config1 = Config::ReadPOSCAR(filename1, false);
  auto config2 = Config::ReadPOSCAR(filename2, false);
  auto index_and_distance_squared_vector = GetSortedIndexDistanceSquaredHelper(config1, config2);
  // we assume the atom that moves the most is the jump atom
  const int jump_atom_index = index_and_distance_squared_vector.back().first;

  Config config(config1);
  config.AppendAtom({config1.GetNumAtoms(), 0, "X",
                     config2.GetAtomList()[jump_atom_index].GetRelativePosition()});
  config.UpdateNeighbors();
  const Atom &jump_atom = config.GetAtomList()[jump_atom_index];
  const Atom &vacancy = config.GetAtomList().back();
  // Check and warning
  if (jump_atom.GetFirstNearestNeighborList().size() != Al_const::kNumFirstNearestNeighbors ||
      vacancy.GetFirstNearestNeighborList().size() != Al_const::kNumFirstNearestNeighbors ||
      jump_atom.GetSecondNearestNeighborList().size() != Al_const::kNumSecondNearestNeighbors ||
      vacancy.GetSecondNearestNeighborList().size() != Al_const::kNumSecondNearestNeighbors ||
      jump_atom.GetThirdNearestNeighborList().size() != Al_const::kNumThirdNearestNeighbors ||
      vacancy.GetThirdNearestNeighborList().size() != Al_const::kNumThirdNearestNeighbors) {
    std::cerr << "# number of neighbors warning\n# ";

    std::cerr << jump_atom.GetFirstNearestNeighborList().size() << ' '
              << vacancy.GetFirstNearestNeighborList().size() << ' '
              << jump_atom.GetSecondNearestNeighborList().size() << ' '
              << vacancy.GetSecondNearestNeighborList().size() << ' '
              << jump_atom.GetThirdNearestNeighborList().size() << ' '
              << vacancy.GetThirdNearestNeighborList().size() << '\n';
    // return;
  }

  std::unordered_set<int> atom1_first_neighbors_set(
      jump_atom.GetFirstNearestNeighborList().begin(),
      jump_atom.GetFirstNearestNeighborList().end());
  atom1_first_neighbors_set.erase(vacancy.GetId());

  const std::unordered_set<int> atom1_second_neighbors_set(
      jump_atom.GetSecondNearestNeighborList().begin(),
      jump_atom.GetSecondNearestNeighborList().end());
  const std::unordered_set<int> atom1_third_neighbors_set(
      jump_atom.GetThirdNearestNeighborList().begin(),
      jump_atom.GetThirdNearestNeighborList().end());

  std::unordered_set<int> atom2_first_neighbors_set(
      vacancy.GetFirstNearestNeighborList().begin(),
      vacancy.GetFirstNearestNeighborList().end());
  atom2_first_neighbors_set.erase(jump_atom_index);

  const std::unordered_set<int> atom2_second_neighbors_set(
      vacancy.GetSecondNearestNeighborList().begin(),
      vacancy.GetSecondNearestNeighborList().end());
  const std::unordered_set<int> atom2_third_neighbors_set(
      vacancy.GetThirdNearestNeighborList().begin(),
      vacancy.GetThirdNearestNeighborList().end());

  std::vector<int> near_neighbors_of_jump_pair;
  std::array<std::unordered_set<int>, 10> sub_encode_sets;
  sub_encode_sets[0] = {jump_atom_index};
  for (const int index:atom1_first_neighbors_set) {
    if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
      sub_encode_sets[1].insert(index);
      continue;
    }
    if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
      sub_encode_sets[2].insert(index);
      continue;
    }
    if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
      sub_encode_sets[3].insert(index);
      continue;
    }
    sub_encode_sets[4].insert(index);
  }
  for (const int index:atom2_first_neighbors_set) {
    if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
      continue;
    }
    if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
      sub_encode_sets[2].insert(index);
      continue;
    }
    if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
      sub_encode_sets[3].insert(index);
      continue;
    }
    sub_encode_sets[4].insert(index);
  }

  for (const int index:atom1_second_neighbors_set) {
    if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
      continue;
    }
    if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
      sub_encode_sets[5].insert(index);
      continue;
    }
    if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
      sub_encode_sets[6].insert(index);
      continue;
    }
    sub_encode_sets[7].insert(index);
  }
  for (const int index:atom2_second_neighbors_set) {
    if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
      continue;
    }
    if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
      continue;
    }
    if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
      sub_encode_sets[6].insert(index);
      continue;
    }
    sub_encode_sets[7].insert(index);
  }

  for (const int index:atom1_third_neighbors_set) {
    if (atom2_first_neighbors_set.find(index) != atom2_first_neighbors_set.end()) {
      continue;
    }
    if (atom2_second_neighbors_set.find(index) != atom2_second_neighbors_set.end()) {
      continue;
    }
    if (atom2_third_neighbors_set.find(index) != atom2_third_neighbors_set.end()) {
      sub_encode_sets[8].insert(index);
      continue;
    }
    sub_encode_sets[9].insert(index);
  }
  for (const int index:atom2_third_neighbors_set) {
    if (atom1_first_neighbors_set.find(index) != atom1_first_neighbors_set.end()) {
      continue;
    }
    if (atom1_second_neighbors_set.find(index) != atom1_second_neighbors_set.end()) {
      continue;
    }
    if (atom1_third_neighbors_set.find(index) != atom1_third_neighbors_set.end()) {
      continue;
    }
    sub_encode_sets[9].insert(index);
  }

  std::unordered_map<int, double>
      index_and_distance_squared_map(index_and_distance_squared_vector.begin(),
                                     index_and_distance_squared_vector.end());
  double sum_of_distance = 0.0;
  for (const auto &sub_encode_set:sub_encode_sets) {
    for (const auto &index:sub_encode_set) {
      sum_of_distance+= index_and_distance_squared_map[index];
    }
    std::cout << std::sqrt(sum_of_distance) << " ";
  }

  std::cout << std::sqrt(std::accumulate(index_and_distance_squared_vector.begin(),
                                         index_and_distance_squared_vector.end(), 0.0,
                                         [](double sum, std::pair<int, double> pair) {
                                           return sum + pair.second;
                                         }));
}
} // namespace kn::NebDistance