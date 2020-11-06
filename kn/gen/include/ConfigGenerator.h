#ifndef KN_KN_NEB_INCLUDE_CONFIGGENERATOR_H_
#define KN_KN_NEB_INCLUDE_CONFIGGENERATOR_H_
#include <filesystem>

#include "Config.h"

namespace gen {
class ConfigGenerator {
  public:
    ConfigGenerator(double lattice_constant,
                    Factor_t factors,
                    std::string solvent_element,
                    std::set<std::string> element_set,
                    const std::filesystem::path &pot_folder_path);
    virtual ~ConfigGenerator();
    virtual void CreateConfigs() const = 0;

  protected:
    void PrepareVASPFiles(const cfg::Config &reference_config,
                          const std::filesystem::path &file_path) const;

    double lattice_constant_;
    Factor_t factors_;
    std::string solvent_element_;
    std::set<std::string> element_set_;
    std::string pot_folder_path_;

    mutable std::mt19937_64 generator_;

};

} // namespace gen
#endif //KN_KN_NEB_INCLUDE_CONFIGGENERATOR_H_
