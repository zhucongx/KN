#include "BarrierPredictor.h"

#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace kmc {
BarrierPredictor::BarrierPredictor(
    const std::string &predictor_filename,
    const cfg::Config &reference_config,
    std::unordered_map<std::string, double> type_category_hashmap)
    : mapping_(ansys::ClusterExpansion::GetAverageClusterParametersMapping(reference_config)),
      type_category_hashmap_(std::move(type_category_hashmap)) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters] : all_parameters.items()) {

    element_parameters_hashmap_[element] = Element_Parameters{
        parameters.at("mu_x"),
        parameters.at("transform_matrix"),
        parameters.at("theta"),
        parameters.at("mean_y")};
  }
}
std::pair<double, double> BarrierPredictor::GetBarrierAndDiff(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) const {

  auto[forward_encode_list, backward_encode_list] =
  ansys::ClusterExpansion::GetAverageClusterParametersForwardAndBackwardFromMap(
      config, jump_pair, type_category_hashmap_, mapping_);
  const auto &element_parameter =
      element_parameters_hashmap_.at(config.GetAtomList().at(jump_pair.second).GetType());
  const auto &mu_x = element_parameter.mu_x;
  const auto &transform_matrix = element_parameter.transform_matrix;
  const auto &theta = element_parameter.theta;
  const auto mean_y = element_parameter.mean_y;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    forward_encode_list[i] -= mu_x[i];
    backward_encode_list[i] -= mu_x[i];
  }
  double forward_barrier = 0;
  double backward_barrier = 0;

  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot_forward = 0;
    double pca_dot_backward = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot_forward += forward_encode_list[i] * transform_matrix[j][i];
      pca_dot_backward += backward_encode_list[i] * transform_matrix[j][i];
    }
    forward_barrier += pca_dot_forward * theta[j];
    backward_barrier += pca_dot_backward * theta[j];
  }
  forward_barrier -= mean_y;
  backward_barrier -= mean_y;
  // std::cerr << config.GetAtomList().at(jump_pair.second).GetType() << "  ";

  // std::cerr << forward_barrier << "\n";
  // static std::mt19937_64 generator(static_cast<unsigned long long int>(
  //                                      std::chrono::system_clock::now().time_since_epoch().count()));
  // static std::uniform_real_distribution<double> distribution(0, 1e-1);
  // auto non_neg_forward = std::max(forward_barrier, 1e-3);
  return {forward_barrier,
          forward_barrier - backward_barrier};
}
} // namespace kmc