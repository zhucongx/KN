#ifndef KN_KN_NEB_INCLUDE_SIZEMISFITGENERATOR_H_
#define KN_KN_NEB_INCLUDE_SIZEMISFITGENERATOR_H_
#include "ConfigGenerator.h"

namespace neb {
class SizeMisfitGenerator : public ConfigGenerator {
  public:
    SizeMisfitGenerator(double lattice_constant,
                        const Factor_t &factors,
                        const std::string &solvent_element,
                        const std::set<std::string> &element_list,
                        const std::filesystem::path &pot_folder_path);
    virtual ~SizeMisfitGenerator();
    void CreateConfigs() const override;
  private:
};

} // namespace neb
#endif //KN_KN_NEB_INCLUDE_SIZEMISFITGENERATOR_H_
