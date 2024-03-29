#include "ClusterConfigGenerator.h"

#include <utility>
#include "ClusterExpansion.h"
namespace gen {
ClusterConfigGenerator::ClusterConfigGenerator(double lattice_constant,
                                               const Factor_t &factors,
                                               const std::string &solvent_element,
                                               const std::set<std::string> &element_set,
                                               const std::filesystem::path &pot_folder_path)
    : ConfigGenerator(lattice_constant, factors, solvent_element, element_set, pot_folder_path) {}
ClusterConfigGenerator::~ClusterConfigGenerator() = default;
static std::vector<size_t> GetEquivalentSingletIndexVector(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) {
  // Get first, second, third nearest neighbors of the jump pairs
  std::unordered_set<size_t>
      atom_id_hashset = GetFirstAndSecondThirdNeighborsSetOfJumpPair(config, atom_id_jump_pair);
  atom_id_hashset.erase(atom_id_jump_pair.first);
  atom_id_hashset.erase(atom_id_jump_pair.second);

  std::vector<cfg::Atom> atom_list;
  atom_list.reserve(atom_id_hashset.size());

  const auto move_distance = Vector_t{0.5, 0.5, 0.5} - GetPairCenter(config, atom_id_jump_pair);
  for (auto id : atom_id_hashset) {
    cfg::Atom atom = config.GetAtomList()[id];

    // move to center
    auto relative_position = atom.GetRelativePosition();
    relative_position += move_distance;
    relative_position -= ElementFloor(relative_position);
    atom.SetRelativePosition(relative_position);
    atom_list.push_back(std::move(atom));
  }
  RotateAtomVector(atom_list, GetPairRotationMatrix(config, atom_id_jump_pair));

  auto IsSmallerSymmetrically = [](const auto &lhs, const auto &rhs) -> bool {
    const auto &relative_position_lhs = lhs.GetRelativePosition();
    const auto &relative_position_rhs = rhs.GetRelativePosition();

    const double diff_x = relative_position_lhs[kXDimension] - relative_position_rhs[kXDimension];
    if (diff_x < -kEpsilon)
      return true;
    if (diff_x > kEpsilon)
      return false;

    const double diff_y = std::abs(relative_position_lhs[kYDimension] - 0.5)
        - std::abs(relative_position_rhs[kYDimension] - 0.5);
    if (diff_y < -kEpsilon)
      return true;
    if (diff_y > kEpsilon)
      return false;

    return (std::abs(relative_position_lhs[kZDimension] - 0.5)
        < std::abs(relative_position_rhs[kZDimension] - 0.5) - kEpsilon);
  };
  std::set<cfg::Atom, decltype(IsSmallerSymmetrically)> atom_set;
  for (auto &atom:atom_list) {
    atom_set.insert(std::move(atom));
  }

  std::vector<size_t> result;
  result.reserve(atom_set.size());
  for (const auto &atom:atom_set) {
    result.emplace_back(atom.GetId());
  }
  return result;
}

void ClusterConfigGenerator::CreateSingletsConfigs() const {
  auto base_config = cfg::GenerateFCC(lattice_constant_, solvent_element_, factors_);
  base_config.ChangeAtomTypeAt(0, "X");
  // For convenience choose jump pair 0 and 1, where X at 0 and jump_atom at 1
  const std::pair<size_t, size_t> atom_id_jump_pair{0, 1};
  auto singlet_id_vector = GetEquivalentSingletIndexVector(base_config, atom_id_jump_pair);

  std::ofstream ofs("log.txt", std::ofstream::out);
  size_t count = 0;
  for (const auto &jump_type : element_set_) {
    auto reference_config = base_config;
    reference_config.ChangeAtomTypeAt(1, jump_type);
    for (const auto &singlet_type:element_set_) {
      if (singlet_type == solvent_element_)
        continue;
      for (auto singlet_id:singlet_id_vector) {
        auto config_start = reference_config;
        config_start.ChangeAtomTypeAt(singlet_id, singlet_type);
        auto config_end = config_start;
        cfg::AtomsJump(config_end, atom_id_jump_pair);

        cfg::Config::WriteConfig(config_start, std::to_string(count) + ".cfg", 0);

        ofs << "config " << count << " end 0 pair: "
            << atom_id_jump_pair.first << ' ' << atom_id_jump_pair.second << '\n';
        std::filesystem::path config_path("config" + std::to_string(count));
        // Generate start files
        std::filesystem::path start_path(config_path / "s");
        std::filesystem::create_directories(start_path);
        cfg::Config::WriteConfig(config_start, start_path / "start.cfg", 0);
        config_start.Perturb(generator_);
        cfg::Config::WritePOSCAR(config_start, start_path / "POSCAR", false);
        PrepareVASPFiles(config_start, start_path);
        // Generate end files
        std::filesystem::path end_path(config_path / "e_0");
        std::filesystem::create_directories(end_path);
        cfg::Config::WriteConfig(config_end, end_path / "end.cfg", 0);
        config_end.Perturb(generator_);
        cfg::Config::WritePOSCAR(config_end, end_path / "POSCAR", false);
        PrepareVASPFiles(config_end, end_path);
        ++count;
      }
    }
  }
}
void ClusterConfigGenerator::CreateConfigs() const {
// Find equivalent atoms, pairs, triplets
  CreateSingletsConfigs();
}

} // namespace gen::cluster_model