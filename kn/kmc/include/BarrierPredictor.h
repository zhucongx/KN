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
                     const std::set<std::string> &type_set);
    virtual ~BarrierPredictor();

    [[nodiscard]] virtual std::pair<double, double> GetBarrierAndDiff(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &jump_pair) const;
  protected:
    [[nodiscard]] double GetBarrierFromEncode(
        const std::string &element_type,
        const std::vector<std::string>& encode_list) const;
  private:
    std::unordered_map<std::string, Element_Parameters> element_parameters_hashmap_{};
    const std::vector<std::vector<std::vector<size_t>>> mapping_;
    const std::unordered_map<std::string, std::vector<double>> one_hot_encode_hash_map_;
    const size_t num_of_elements_;

};
} // namespace kmc
#endif //KN_KN_KMC_INCLUDE_BARRIERPREDICTOR_H_
