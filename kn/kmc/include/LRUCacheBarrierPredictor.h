#ifndef KN_KN_KMC_INCLUDE_LRUCACHEBARRIERPREDICTOR_H_
#define KN_KN_KMC_INCLUDE_LRUCACHEBARRIERPREDICTOR_H_
#include "BarrierPredictor.h"
#include <boost/functional/hash.hpp>

namespace kmc {
class LRUCacheBarrierPredictor : BarrierPredictor {
  public:
    LRUCacheBarrierPredictor(const std::string &predictor_filename,
                             const cfg::Config &reference_config,
                             const std::set<std::string> &type_set,
                             size_t cache_size);
    ~LRUCacheBarrierPredictor() override;

    [[nodiscard]] std::pair<double, double> GetBarrierAndDiff(
        const cfg::Config &config,
        const std::pair<size_t, size_t> &jump_pair) const override;
  private:
    void Add(const std::vector<std::string> &key, double value) const;
    size_t cache_size_;
    mutable std::list<std::pair<std::vector<std::string>, double>> cache_list_{};
    mutable std::unordered_map<std::vector<std::string>,
                               std::list<std::pair<std::vector<std::string>, double>>::iterator,
                               boost::hash<std::vector<std::string>>> hashmap_{};
#ifndef NDEBUG
  public:
    mutable size_t count_{0};
#endif
};
} // namespace kmc


#endif //KN_KN_KMC_INCLUDE_LRUCACHEBARRIERPREDICTOR_H_
