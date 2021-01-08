#ifndef KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
#define KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
#include "ClusterExpansion.h"


namespace kmc {
struct Element_Parameters {
  std::vector<double> mu_x;
  std::vector<std::vector<double>> transform_matrix;
  std::vector<double> theta;
  double mean_y;
};

class BarrierPredictor {
  public:
    BarrierPredictor(const std::string &predictor_filename,
                     const cfg::Config &reference_config,
                     std::unordered_map<std::string, double> type_category_hashmap);

    [[nodiscard]] std::pair<double, double> GetBarrierAndDiff(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &jump_pair) const;
  private:
    std::unordered_map<std::string, Element_Parameters> element_parameters_hashmap_{};
    const ansys::ClusterExpansion::Cluster_Map_t mapping_;
    const std::unordered_map<std::string, double> type_category_hashmap_;
};
} // namespace kmc
#endif //KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
