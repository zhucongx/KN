#ifndef KN_KN_KMC_INCLUDE_BEPBARRIERPREDICTOR_H_
#define KN_KN_KMC_INCLUDE_BEPBARRIERPREDICTOR_H_
#include "BarrierPredictor.h"
namespace kmc {
class BEPBarrierPredictor : BarrierPredictor {
  public:
    BEPBarrierPredictor(const std::string &predictor_filename,
                        const cfg::Config &reference_config,
                        const std::set<std::string> &type_set);
    ~BEPBarrierPredictor() override;

    [[nodiscard]] std::pair<double, double> GetBarrierAndDiff(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &atom_id_jump_pair) const override;
  private:
    const std::map<std::string, double> standard_e0_{{"Al", 0.52}, {"Mg", 0.45}, {"Zn", 0.34}};
};
} // namespace kmc

#endif //KN_KN_KMC_INCLUDE_BEPBARRIERPREDICTOR_H_
