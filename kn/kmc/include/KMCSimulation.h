#ifndef KN_KN_KMC_INCLUDE_KMCSIMULATION_H_
#define KN_KN_KMC_INCLUDE_KMCSIMULATION_H_

#include "KMCEvent.h"
#include "LRUCacheBarrierPredictor.h"

namespace kmc {
class KMCSimulation {
  public:
    KMCSimulation(cfg::Config config,
                  unsigned long long int log_dump_steps,
                  unsigned long long int config_dump_steps,
                  unsigned long long int maximum_number,
                  const std::set<std::string> &type_set,
                  unsigned long long int steps,
                  double energy,
                  double time,
                  const std::string &json_parameters_filename,
                  size_t lru_size);
    virtual ~KMCSimulation();
    virtual void CheckTimeAndFix(double one_step_time);
    virtual void Simulate();
  protected:
    void BuildEventListSerial();
    void BuildEventListParallel();
    size_t SelectEvent() const;

    // simulation parameters
    cfg::Config config_;
    const unsigned long long int log_dump_steps_;
    const unsigned long long int config_dump_steps_;
    const unsigned long long int maximum_number_;

    // simulation statistics
    unsigned long long int steps_;
    double energy_;
    double time_;
    const size_t vacancy_index_;
    double one_step_barrier_{0.0};
    double one_step_change_{0.0};

    // helpful properties
    double total_rate_{0.0};
    size_t mpi_rank_;
    size_t mpi_size_;

    std::vector<KMCEvent> event_list_{};
    const LRUCacheBarrierPredictor lru_cache_barrier_predictor_;
    mutable std::mt19937_64 generator_;
};
} // namespace kmc


#endif //KN_KN_KMC_INCLUDE_KMCSIMULATION_H_
