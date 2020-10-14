#ifndef KN_KN_NEB_INCLUDE_CONFIGGENERATOR_H_
#define KN_KN_NEB_INCLUDE_CONFIGGENERATOR_H_
#include "Config.h"
namespace neb {
class ConfigGenerator {
  public:
    ConfigGenerator(double lattice_constant,
                    Factor_t factors,
                    std::filesystem::path solvent_element,
                    std::set<std::string> element_list,
                    std::filesystem::path pot_folder_path);
    virtual void CreateConfigs() const = 0;

  protected:
    void PrepareVASPFiles(const cfg::Config &reference_config,
                          const std::filesystem::path &file_path) const;

    double lattice_constant_;
    Factor_t factors_;
    std::string solvent_element;
    std::set<std::string> element_set_;
    std::string pot_folder_path_;

    mutable std::mt19937_64 generator_;

};

} // namespace cfg
#endif //KN_KN_NEB_INCLUDE_CONFIGGENERATOR_H_
