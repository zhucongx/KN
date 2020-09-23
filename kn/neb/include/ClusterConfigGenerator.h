#ifndef KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
#define KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
#include "ConfigGenerator.h"
namespace neb {
// This class generate some singular atom, pairs triplets as criterion
class ClusterConfigGenerator : public ConfigGenerator {
  public:
    ClusterConfigGenerator(double lattice_constant,
                           const Factor_t &factors,
                           const std::string &solvent_element,
                           const std::set<std::string> &element_list,
                           const std::string &pot_folder_path);
    void CreateConfigs() override;
  private:
    void CreateSingletsConfigs();

};

} // namespace neb
#endif //KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
