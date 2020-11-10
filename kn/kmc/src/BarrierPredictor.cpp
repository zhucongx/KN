#include "BarrierPredictor.h"
#include <boost/range/combine.hpp>

#include <utility>
namespace kmc {
BarrierPredictor::BarrierPredictor(
    const cfg::Config &reference_config,
    std::unordered_map<std::string, double> type_category_hashmap)
    : type_category_hashmap_(std::move(type_category_hashmap)) {
  mapping_ = ansys::ClusterExpansion::GetAverageClusterParametersMapping(reference_config);

  auto element_set = reference_config.GetTypeSet();
  for (const auto &element : element_set) {
    std::ifstream ifs(element + ".weight", std::ifstream::in);
    if (!ifs.is_open()) {
      std::cout << "Cannot open " << element << ".weight\n";
      break;
    }
    double weight;
    std::vector<double> weight_vector;
    while (ifs >> weight) {
      weight_vector.emplace_back(weight);
    }
    element_weight_hashmap_.at(element) = weight_vector;
  }
}
std::pair<double, double> BarrierPredictor::GetBarrierAndDiff(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) const {

  auto[forward_encode_list, backward_encode_list] =
  ansys::ClusterExpansion::GetAverageClusterParametersForwardAndBackwardFromMap(
      config, jump_pair, type_category_hashmap_, mapping_);
  const auto &weights
      = element_weight_hashmap_.at(config.GetAtomList()[jump_pair.second].GetType());

  double forward_barrier = 0, backward_barrier = 0;
  for (const auto[forward_encode, backward_encode, weight]
      : boost::combine(forward_encode_list, backward_encode_list, weights)) {
    forward_barrier += forward_encode * weight;
    backward_barrier += backward_encode * weight;
  }

  return {forward_barrier, forward_barrier - backward_barrier};
}
} // namespace kmc