#ifndef KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
#define KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
#include "ConfigGenerator.h"
namespace neb {
// This class generate some singular atom, pairs triplets as criterion
class ClusterConfigGenerator : public ConfigGenerator {
  public:
    ClusterConfigGenerator(double lattice_constant,
                           const Factor_t &factors,
                           const std::filesystem::path &solvent_element,
                           const std::set<std::string> &element_list,
                           const std::filesystem::path &pot_folder_path);
    void CreateConfigs() const override;
  private:
    void CreateSingletsConfigs() const;

};

} // namespace neb
#endif //KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
