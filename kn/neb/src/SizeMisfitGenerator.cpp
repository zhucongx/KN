#include "SizeMisfitGenerator.h"
namespace neb {

SizeMisfitGenerator::SizeMisfitGenerator(double lattice_constant,
                                         const Factor_t &factors,
                                         const std::string &solvent_element,
                                         const std::set<std::string> &element_list,
                                         const std::filesystem::path &pot_folder_path)
    : ConfigGenerator(lattice_constant, factors, solvent_element, element_list, pot_folder_path) {}
void SizeMisfitGenerator::CreateConfigs() const {
  auto base_config = cfg::GenerateFCC(lattice_constant_, solvent_element_, factors_);
  for (const auto &element_type : element_set_) {
    std::filesystem::path element_path(element_type);
    auto reference_config = base_config;
    reference_config.ChangeAtomTypeAt(0, element_type);

    double scale = 0.94;
    do {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(2) << scale;
      std::filesystem::path scale_path = element_path / ss.str();
      std::filesystem::create_directories(scale_path);

      auto out_config = reference_config;
      out_config.ScaleWith(scale);
      cfg::Config::WriteConfig(reference_config, scale_path / "supercell.cfg", false);
      out_config.Perturb(generator_);
      cfg::Config::WritePOSCAR(out_config, scale_path / "POSCAR", false);
      PrepareVASPFiles(out_config, scale_path);
      scale += 0.02;
    } while (scale < 1.06 + kEpsilon);
  }
}
} // namespace neb