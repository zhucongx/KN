//
// Created by zhucongx on 1/8/21.
//

#include "LRUCacheBarrierPredictor.h"
namespace kmc {
LRUCacheBarrierPredictor::LRUCacheBarrierPredictor(const std::string &predictor_filename,
                                                   const cfg::Config &reference_config,
                                                   const std::set<std::string> &type_set,
                                                   size_t cache_size)
    : BarrierPredictor(predictor_filename,
                       reference_config,
                       type_set),
      cache_size_(cache_size) {}
LRUCacheBarrierPredictor::~LRUCacheBarrierPredictor() = default;

std::pair<double, double> LRUCacheBarrierPredictor::GetBarrierAndDiff(
    const cfg::Config &config,
    const std::pair<size_t, size_t> &jump_pair) const {

  const auto &element_type = config.GetAtomList().at(jump_pair.second).GetType();
  const auto[forward_encode_list, backward_encode_list] =
  ansys::ClusterExpansion::GetForwardAndBackwardEncode(config, jump_pair);
  double forward_e0, backward_e0;
  auto it1 = hashmap_.find(forward_encode_list);
  if (it1 == hashmap_.end()) {
    forward_e0 = GetE0FromEncode(element_type, forward_encode_list);
    Add(forward_encode_list, forward_e0);
  } else {
    cache_list_.splice(cache_list_.begin(), cache_list_, it1->second);
    forward_e0 = it1->second->second;
  }

  auto it2 = hashmap_.find(backward_encode_list);
  if (it2 == hashmap_.end()) {
    backward_e0 = GetE0FromEncode(element_type, backward_encode_list);
    Add(backward_encode_list, backward_e0);
  } else {
#ifndef NDEBUG
    ++count_;
#endif
    cache_list_.splice(cache_list_.begin(), cache_list_, it2->second);
    backward_e0 = it2->second->second;
  }

  auto e0 = 0.5 * (forward_e0 + backward_e0);
  auto dE = GetDEFromConfig(config, jump_pair);
#ifndef NDEBUG
  std::cout << forward_e0 << '\t' << backward_e0 << '\n';
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
void LRUCacheBarrierPredictor::Add(const std::vector<std::string> &key, double value) const {
  auto it = hashmap_.find(key);
  if (it != hashmap_.end()) {
    cache_list_.erase(it->second);
  }
  cache_list_.push_front(std::make_pair(key, value));
  hashmap_[key] = cache_list_.begin();
  if (hashmap_.size() > cache_size_) {
    auto last = cache_list_.rbegin()->first;
    cache_list_.pop_back();
    hashmap_.erase(last);
  }
}

} // namespace kmc
