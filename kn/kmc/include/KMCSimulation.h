#ifndef KN_KN_KMC_INCLUDE_KMCSIMULATION_H_
#define KN_KN_KMC_INCLUDE_KMCSIMULATION_H_

#include <boost/mpi.hpp>

#include "KMCEvent.h"
#include "BarrierPredictor.h"

namespace kmc {
class KMCSimulation {
  public:
    KMCSimulation(cfg::Config config,
                  unsigned long long int log_dump_steps,
                  unsigned long long int config_dump_steps,
                  unsigned long long int maximum_number,
                  std::unordered_map<std::string, double> type_category_hashmap,
                  unsigned long long int steps,
                  double energy,
                  double time);
    // KMCSimulation(cfg::Config config,
    //               unsigned long long int log_dump_steps,
    //               unsigned long long int config_dump_steps,
    //               unsigned long long int maximum_number,
    //               unsigned long long int steps,
    //               double energy,
    //               double time,
    //               std::unordered_map<std::string, double> type_category_hashmap);
    void Simulate();
  protected:
    void BuildEventListSerial();
    void BuildEventListParallel();
    size_t SelectEvent();

    // simulation parameters
    cfg::Config config_;
    const unsigned long long int log_dump_steps;
    const unsigned long long int config_dump_steps;
    const unsigned long long int maximum_number_;

    // statistical information
    unsigned long long int steps_;
    double energy_;
    double time_;
    size_t vacancy_index_;

    // helpful properties
    double total_rate_;
    boost::mpi::environment env_;
    boost::mpi::communicator world_;
    size_t mpi_rank_;
    size_t mpi_size_;

    std::vector<KMCEvent> event_list_;
    BarrierPredictor barrier_predictor_;
    mutable std::mt19937_64 generator_;
};
} // namespace kmc


#endif //KN_KN_KMC_INCLUDE_KMCSIMULATION_H_
