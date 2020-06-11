#ifndef KN_INCLUDE_CONFIGGENERATOR_H_
#define KN_INCLUDE_CONFIGGENERATOR_H_
#include <string>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include "Config.h"

namespace kn {
class ConfigGenerator {
  public:
    using Factor_t = std::array<int, kDimension>;
    ConfigGenerator(double lattice_const,
                    const Factor_t &factors,
                    const std::string &solvent_element,
                    const std::map<std::string, int> &element_count_map,
                    const std::string &pot_folder_path);

    void CreateRandom(int num_configs);
    void CreateSpecific();
  private:
    Config ShuffleConfig(const Config& config) const;
    void PrepareVASPFiles(const std::string &file_path);
    static Config GenerateFCC(double lattice_constant_a,
                              const std::string &element,
                              const Factor_t &factors);
    // static Config GenerateBCC(double lattice_constant_a,
    //                           const std::string &element,
    //                           const Factor_t &factors);
    // static Config GenerateHCP(double lattice_constant_a,
    //                           double lattice_constant_c,
    //                           const std::string &element,
    //                           const Factor_t &factors);
    // Anti-phase Configuration
    static Config GenerateL10(double lattice_constant_a,
                              const std::vector<std::string> &element_list,
                              const Factor_t &factors);
    static Config GenerateL12(double lattice_constant_a,
                              const std::vector<std::string> &element_list,
                              const Factor_t &factors);
    static Config GenerateL10star(double lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const Factor_t &factors);
    static Config GenerateL12star(double lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const Factor_t &factors);
    static Config GenerateZ1(double lattice_constant_a,
                             const std::vector<std::string> &element_list,
                             const Factor_t &factors);

  private:
    double lattice_const_;
    Factor_t factors_;
    std::string solvent_element_;
    // element name -- the number of this element in the config
    std::map<std::string, int> element_count_map_;
    // vector of elements
    std::vector<std::string> element_list_;
    // vector of atoms' types
    std::vector<std::string> atom_index_list_;
    std::string pot_folder_path_;

    mutable std::mt19937_64 generator_;
};
} // namespace kn

#endif //KN_INCLUDE_CONFIGGENERATOR_H_
