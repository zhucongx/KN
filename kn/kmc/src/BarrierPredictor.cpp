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
    : initialized_hashmap_(ansys::BondCounting::InitializeHashMap(type_set)),
      mapping_(ansys::ClusterExpansion::GetAverageClusterParametersMapping(reference_config)),
      one_hot_encode_hash_map_(ansys::ClusterExpansion::GetOneHotEncodeHashmap(type_set)),
      num_of_elements_(type_set.size()) {
  std::ifstream ifs(predictor_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;
  for (const auto &[element, parameters]: all_parameters.items()) {
    if (element == "Bond") {
      theta_ = std::vector<double>(parameters.at("theta"));
      continue;
    }
    element_parameters_hashmap_[element] = Element_Parameters{
        parameters.at("mu_x"),
        parameters.at("sigma_x"),
        parameters.at("mu_y"),
        parameters.at("sigma_y"),
        parameters.at("U"),
        parameters.at("theta"),
    };
  }
}
BarrierPredictor::~BarrierPredictor() = default;
double BarrierPredictor::GetE0FromEncode(
    const std::string &element_type,
    const std::vector<std::string> &encode_list) const {

  auto one_hot_parameters = ansys::ClusterExpansion::GetOneHotParametersFromMap(
      encode_list, one_hot_encode_hash_map_, num_of_elements_, mapping_);
  // std::cout << one_hot_parameters.size() << std::endl;

  const auto &element_parameters = element_parameters_hashmap_.at(element_type);
  const auto &mu_x = element_parameters.mu_x;
  const auto &sigma_x = element_parameters.sigma_x;
  const auto mu_y = element_parameters.mu_y;
  const auto sigma_y = element_parameters.sigma_y;

  const auto &U = element_parameters.U;
  const auto &theta = element_parameters.theta;

  const size_t old_size = mu_x.size();
  const size_t new_size = theta.size();

  for (size_t i = 0; i < old_size; ++i) {
    one_hot_parameters[i] -= mu_x[i];
    one_hot_parameters[i] /= sigma_x[i];
  }
  double e0 = 0;
  for (size_t j = 0; j < new_size; ++j) {
    double pca_dot = 0;
    for (size_t i = 0; i < old_size; ++i) {
      pca_dot += one_hot_parameters[i] * U[j][i];
    }
    e0 += pca_dot * theta[j];
  }
  e0 *= sigma_y;
  e0 += mu_y;
  return e0;
}
double BarrierPredictor::GetDEFromConfig(const cfg::Config &config,
                                         const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  auto bond_change_vector =
      ansys::BondCounting::GetBondChange(config, atom_id_jump_pair, initialized_hashmap_);
  const size_t bond_size = theta_.size();
  double dE = 0;
  for (size_t i = 0; i < bond_size; ++i) {
    dE += theta_[i] * bond_change_vector[i];
  }
  return dE;
  // return 0.0;
}

std::pair<double, double> BarrierPredictor::GetBarrierAndDiff(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &atom_id_jump_pair) const {
  const auto &element_type = config.GetAtomList().at(atom_id_jump_pair.second).GetType();
  auto[forward_encode_list, backward_encode_list] =
  ansys::ClusterExpansion::GetForwardAndBackwardEncode(config, atom_id_jump_pair);

  auto forward_e0 = GetE0FromEncode(element_type, forward_encode_list);
  auto backward_e0 = GetE0FromEncode(element_type, backward_encode_list);
  auto e0 = 0.5 * (forward_e0 + backward_e0);

  auto dE = GetDEFromConfig(config, atom_id_jump_pair);

  // std::cerr << config.GetAtomList().at(atom_id_jump_pair.second).GetType() << "  ";

  // std::cerr << forward_barrier << "\n";
  // static std::mt19937_64 generator(static_cast<unsigned long long int>(
  //                                      std::chrono::system_clock::now().time_since_epoch().count()));
  // static std::uniform_real_distribution<double> distribution(0, 1e-1);
  // auto non_neg_forward = std::max(forward_barrier, 1e-3);
  return {e0 + dE / 2, dE};
}
} // namespace kmc