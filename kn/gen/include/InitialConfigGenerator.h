#ifndef KN_KN_GEN_INCLUDE_INITIALCONFIGGENERATOR_H_
#define KN_KN_GEN_INCLUDE_INITIALCONFIGGENERATOR_H_
#include <string>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include "Config.h"
namespace gen {
class InitialConfigGenerator {
  public:
    InitialConfigGenerator(double lattice_const,
                    const Factor_t &factors,
                    std::string solvent_element,
                    std::map<std::string, size_t> element_count_map);
    static cfg::Config EmbedToLarge(double lattice_const,
                             const Factor_t &new_factors,
                             const Factor_t &old_factors,
                             const cfg::Config &small_config,
                             std::map<std::string, size_t> element_number_map);
    cfg::Config ShuffleConfig(const cfg::Config& config) const;

    // Anti-phase Configuration
    static cfg::Config GenerateL10(double lattice_constant_a,
                                   const std::vector<std::string> &element_list,
                                   const Factor_t &factors);
    static cfg::Config GenerateL12(double lattice_constant_a,
                                   const std::vector<std::string> &element_list,
                                   const Factor_t &factors);
    static cfg::Config GenerateL10star(double lattice_constant_a,
                                       const std::vector<std::string> &element_list,
                                       const Factor_t &factors);
    static cfg::Config GenerateL12star(double lattice_constant_a,
                                       const std::vector<std::string> &element_list,
                                       const Factor_t &factors);
    static cfg::Config GenerateZ1(double lattice_constant_a,
                                  const std::vector<std::string> &element_list,
                                  const Factor_t &factors);
  private:
    double lattice_const_;
    Factor_t factors_;
    std::string solvent_element_;
    // element name -- the number of this element in the config
    std::map<std::string, size_t> element_count_map_;
    // vector of elements
    std::vector<std::string> element_list_;
    // vector of atoms' types
    std::vector<std::string> atom_index_list_;

    mutable std::mt19937_64 generator_;
};
} // namespace cfg
#endif //KN_KN_GEN_INCLUDE_INITIALCONFIGGENERATOR_H_
