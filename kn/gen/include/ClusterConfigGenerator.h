#ifndef KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
#define KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
#include "ConfigGenerator.h"
namespace gen {
// This class generate some singular atom, pairs triplets as criterion
class ClusterConfigGenerator : public ConfigGenerator {
  public:
    ClusterConfigGenerator(double lattice_constant,
                           const Factor_t &factors,
                           const std::string &solvent_element,
                           const std::set<std::string> &element_set,
                           const std::filesystem::path &pot_folder_path);
    ~ClusterConfigGenerator() override;
    void CreateConfigs() const override;
  private:
    void CreateSingletsConfigs() const;

};

} // namespace gen
#endif //KN_KN_NEB_INCLUDE_CLUSTERCONFIGGENERATOR_H_
