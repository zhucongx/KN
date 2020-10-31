#include "SizeMisfitGenerator.h"
namespace neb {
SizeMisfitGenerator::SizeMisfitGenerator(double lattice_constant,
                                         const Factor_t &factors,
                                         const std::string &solvent_element,
                                         const std::set<std::string> &element_list,
                                         const std::filesystem::path &pot_folder_path)
    : ConfigGenerator(lattice_constant, factors, solvent_element, element_list, pot_folder_path) {}
SizeMisfitGenerator::~SizeMisfitGenerator() = default;

static void OverwriteINCAR(const std::filesystem::path &path) {
  std::ofstream ofs(path / "INCAR", std::ofstream::out);
  ofs << "NWRITE = 2       \n"
      << "GGA = PE         \n"
      << "\n"
      << "PREC   = Accurate\n"
      << "ISYM   = 2       \n"
      << "SYMPREC = 1e-8   \n"
      << "NELM   = 240     \n"
      << "NELMIN = 4       \n"
      << "\n"
      << "NSW    = 10000   \n"
      << "IBRION = 2       \n"
      << "POTIM  = 0.5     \n"
      << "ISIF   = 2       \n"
      << "\n"
      << "ISMEAR = 1       \n"
      << "SIGMA  = 0.05    \n"
      << "\n"
      << "LREAL  = AUTO    \n"
      << "ENCUT  = 750.00  \n"
      << "ENAUG  = 800.00  \n"
      << "EDIFF  = 1e-7    \n"
      << "ISPIN  = 1       \n"
      << "\n"
      << "LWAVE  = .FALSE. \n"
      << "LCHARG = .TRUE.  \n"
      << "                 \n"
      << "NPAR   = 4       \n";
}

void SizeMisfitGenerator::CreateConfigs() const {
  auto base_config = cfg::GenerateFCC(lattice_constant_, solvent_element_, factors_);
  for (const auto &element_type : element_set_) {
    std::filesystem::path element_path(element_type);
    auto reference_config = base_config;
    reference_config.ChangeAtomTypeAt(0, element_type);

    for (const double scale : {0.940, 0.950, 0.960, 0.970, 0.980, 0.990, 0.992, 0.994, 0.996, 0.998,
                               1.000, 1.002, 1.004, 1.006, 1.008, 1.010, 1.020, 1.030, 1.040, 1.050,
                               1.060}) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(3) << scale;
      std::filesystem::path scale_path = element_path / ss.str();
      std::filesystem::create_directories(scale_path);

      auto out_config = reference_config;
      out_config.ScaleWith(scale);
      cfg::Config::WriteConfig(reference_config, scale_path / "supercell.cfg", false);
      // out_config.Perturb(generator_);
      cfg::Config::WritePOSCAR(out_config, scale_path / "POSCAR", false);
      PrepareVASPFiles(out_config, scale_path);
      OverwriteINCAR(scale_path);
    }
  }
}

} // namespace neb