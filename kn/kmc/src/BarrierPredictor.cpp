#include "BarrierPredictor.h"

#include <utility>
#include <boost/range/combine.hpp>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace kmc {
BarrierPredictor::BarrierPredictor(
    const std::string &predictor_filename,
    const cfg::Config &reference_config,
    const std::set<std::string> &type_set)
    : mapping_(ansys::ClusterExpansion::GetAverageClusterParametersMapping(reference_config)),
      one_hot_encode_hash_map_(ansys::ClusterExpansion::GetOneHotEncodeHashmap(type_set)),
      num_of_elements_(type_set.size()) {
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
BarrierPredictor::~BarrierPredictor() = default;
double BarrierPredictor::GetBarrierFromEncode(
    const std::string &element_type,
    const std::vector<std::string>& encode_list) const{

  auto one_hot_parameters = ansys::ClusterExpansion::GetOneHotParametersFromMap(
      encode_list, one_hot_encode_hash_map_, num_of_elements_, mapping_);

  const auto &element_parameter = element_parameters_hashmap_.at(element_type);
  const auto &mu_x = element_parameter.mu_x;
  const auto &transform_matrix = element_parameter.transform_matrix;
  const auto &theta = element_parameter.theta;
  const auto mean_y = element_parameter.mean_y;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    one_hot_parameters[i] -= mu_x[i];
  }
  double barrier = 0;

  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += one_hot_parameters[i] * transform_matrix[j][i];
    }
    barrier += pca_dot * theta[j];
  }
  barrier += mean_y;

  return barrier;
}

std::pair<double, double> BarrierPredictor::GetBarrierAndDiff(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) const {
  const auto &element_type = config.GetAtomList().at(jump_pair.second).GetType();
  auto[forward_encode_list, backward_encode_list] =
  ansys::ClusterExpansion::GetForwardAndBackwardEncode(config, jump_pair);

  auto forward_barrier = GetBarrierFromEncode(element_type, forward_encode_list);
  auto backward_barrier = GetBarrierFromEncode(element_type, backward_encode_list);

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