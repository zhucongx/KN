#include "BEPBarrierPredictor.h"

namespace kmc {

kmc::BEPBarrierPredictor::BEPBarrierPredictor(const std::string &predictor_filename,
                                              const cfg::Config &reference_config,
                                              const std::set<std::string> &type_set)
    : BarrierPredictor(predictor_filename, reference_config, type_set) {}
kmc::BEPBarrierPredictor::~BEPBarrierPredictor() {

}
std::pair<double, double> BEPBarrierPredictor::GetBarrierAndDiff(const cfg::Config &config,
                                                                 const std::pair<size_t,
                                                                                 size_t> &jump_pair) const {
  const auto &element_type = config.GetAtomList().at(jump_pair.second).GetType();
  auto e0 = standard_e0_[element_type];

  auto dE = GetDEFromConfig(config, jump_pair);
#ifndef NDEBUG
  std::cout << dE << '\n';
#endif
  // std::cerr << config.GetAtomList().at(jump_pair.second).GetType() << "  ";

  // std::cerr << forward_barrier << "\n";
  // static std::mt19937_64 generator(static_cast<unsigned long long int>(
  //                                      std::chrono::system_clock::now().time_since_epoch().count()));
  // static std::uniform_real_distribution<double> distribution(0, 1e-1);
  // auto non_neg_forward = std::max(forward_barrier, 1e-3);
  return {e0 + dE / 2, dE};
}
} // namespace kmc