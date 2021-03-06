#ifndef KN_KN_NEB_INCLUDE_SIZEMISFITGENERATOR_H_
#define KN_KN_NEB_INCLUDE_SIZEMISFITGENERATOR_H_
#include "ConfigGenerator.h"

namespace gen {
class SizeMisfitGenerator : public ConfigGenerator {
  public:
    SizeMisfitGenerator(double lattice_constant,
                        const Factor_t &factors,
                        const std::string &solvent_element,
                        const std::set<std::string> &element_list,
                        const std::filesystem::path &pot_folder_path);
    ~SizeMisfitGenerator() override;
    void CreateConfigs() const override;
  private:
};

} // namespace gen
#endif //KN_KN_NEB_INCLUDE_SIZEMISFITGENERATOR_H_
