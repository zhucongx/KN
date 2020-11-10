#ifndef KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
#define KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
#include "ClusterExpansion.h"

namespace kmc {
class BarrierPredictor {
  public:
    BarrierPredictor(const cfg::Config &reference_config,
                     std::unordered_map<std::string, double> type_category_hashmap);

    [[nodiscard]] std::pair<double, double> GetBarrierAndDiff(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &jump_pair) const;
  private:
    std::unordered_map<std::string, std::vector<double>> element_weight_hashmap_;
    ansys::ClusterExpansion::Cluster_Map_t mapping_;
    const std::unordered_map<std::string, double> type_category_hashmap_;
};
} // namespace kmc
#endif //KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
