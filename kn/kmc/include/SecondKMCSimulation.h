#ifndef KN_KN_KMC_INCLUDE_SECONDKMCSIMULATION_H_
#define KN_KN_KMC_INCLUDE_SECONDKMCSIMULATION_H_

#include "KMCEvent.h"
#include "LRUCacheBarrierPredictor.h"
#include <mpi.h>
/*This method may be wrong. Do not use*/
namespace kmc {
class SecondKMCSimulation {
  public:
    SecondKMCSimulation(cfg::Config config,
                        unsigned long long int log_dump_steps,
                        unsigned long long int config_dump_steps,
                        unsigned long long int maximum_number,
                        const std::set<std::string> &type_set,
                        unsigned long long int steps,
                        double energy,
                        double time,
                        const std::string &json_parameters_filename,
                        size_t lru_size);
    virtual ~SecondKMCSimulation();
    void Simulate();

  protected:
    void BuildEventListParallel();
    std::pair<size_t, size_t> BuildProbabilityListParallel(double &first_probability,
                                                           double &first_energy_change);
    std::vector<size_t> GetSecondNeighborsIndexes();
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
    double first_total_rate_{0.0}, second_total_rate_{0.0};
    int world_rank_{-1}, first_group_rank_{-1}, second_group_rank_{-1};

    MPI_Group world_group_, first_group_, second_group_;
    MPI_Comm first_comm_, second_comm_;

    std::vector<KMCEvent> event_list_{};
    const LRUCacheBarrierPredictor lru_cache_barrier_predictor_;
    mutable std::mt19937_64 generator_;
};

} // namespace kmc
#endif //KN_KN_KMC_INCLUDE_SECONDKMCSIMULATION_H_
