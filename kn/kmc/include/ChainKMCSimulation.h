#ifndef KN_KN_KMC_INCLUDE_CHAINKMCSIMULATION_H_
#define KN_KN_KMC_INCLUDE_CHAINKMCSIMULATION_H_

#include "KMCEvent.h"
#include "LRUCacheBarrierPredictor.h"
#include <mpi.h>
namespace kmc {
//  j -> k -> i ->l
//       |
// current position
class ChainKMCSimulation {
  public:
    ChainKMCSimulation(cfg::Config config,
                       unsigned long long int log_dump_steps,
                       unsigned long long int config_dump_steps,
                       unsigned long long int maximum_number,
                       const std::set<std::string> &type_set,
                       unsigned long long int steps,
                       double energy,
                       double time,
                       const std::string &json_parameters_filename,
                       size_t lru_size);
    virtual ~ChainKMCSimulation();
    virtual void Simulate();

  protected:
    virtual void CheckAndSolveEquilibrium(std::ofstream &ofs){};
    inline void Dump(std::ofstream &ofs);
    KMCEvent GetEventI();
    [[nodiscard]] double BuildEventListParallel();

    std::vector<size_t> GetLIndexes();
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
    double one_step_energy_change_{0.0};
    double one_step_time_change_{0.0};

    // helpful properties
    double total_rate_k_{0.0}, total_rate_i_{0.0};
    int world_rank_{-1}, first_group_rank_{-1}, second_group_rank_{-1};
    // double rij{0}, pij{0};
    size_t previous_j;

    MPI_Group world_group_, first_group_, second_group_;
    MPI_Comm first_comm_, second_comm_;

    std::vector<KMCEvent> event_list_{};
    const LRUCacheBarrierPredictor lru_cache_barrier_predictor_;
    mutable std::mt19937_64 generator_;
};

} // namespace kmc
#endif //KN_KN_KMC_INCLUDE_CHAINKMC_H_
