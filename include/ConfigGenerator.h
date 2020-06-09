#ifndef KN_INCLUDE_CONFIGGENERATOR_H_
#define KN_INCLUDE_CONFIGGENERATOR_H_
#include "Config.h"

namespace kn {
class ConfigGenerator {
  public:
    static Config GenerateFCC(double lattice_constant_a,
                              const std::string &element,
                              const std::array<int, kDimension> &factors);
    static Config GenerateBCC(double lattice_constant_a,
                              const std::string &element,
                              const std::array<int, kDimension> &factors);
    static Config GenerateHCP(double lattice_constant_a,
                              double lattice_constant_c,
                              const std::string &element,
                              const std::array<int, kDimension> &factors);
    // Anti-phase Configuration
    static Config GenerateL10(double lattice_constant_a,
                              const std::vector<std::string> &element_list,
                              const std::array<int, kDimension> &factors);
    static Config GenerateL12(double lattice_constant_a,
                              const std::vector<std::string> &element_list,
                              const std::array<int, kDimension> &factors);
    static Config GenerateL10star(double lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const std::array<int, kDimension> &factors);
    static Config GenerateL12star(double lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const std::array<int, kDimension> &factors);
    static Config GenerateZ1(double lattice_constant_a,
                             const std::vector<std::string> &element_list,
                             const std::array<int, kDimension> &factors);
  private:
    static Config GenerateUnitCell(
        const Matrix33 &basis_matrix,
        const std::vector<std::pair<std::string, Vector3>> &type_position_list);
    static Config Duplicate(const Config &in_config,
                            const std::array<int, kDimension> &factors);
};
} // namespace kn

#endif //KN_INCLUDE_CONFIGGENERATOR_H_