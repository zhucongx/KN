#ifndef KN_INCLUDE_KMCSIMULATION_H_
#define KN_INCLUDE_KMCSIMULATION_H_
#include <fstream>
#include <random>
#include <vector>
#include <boost/mpi.hpp>

#include "Config.h"
#include "KMCEvent.h"
namespace kn {
class KMCSimulation {
  public:
    KMCSimulation(const std::string& cfg_filename,
                  double first_nearest_neighbors_distance = Al_const::kFirstNearestNeighborCutoff,
                  double second_nearest_neighbors_distance = Al_const::kSecondNearestNeighborsCutoff);
    void RunSimulation();
    std::vector<double> CalculateBarrierAndEnergyDifference(const std::pair<int,
                                                                            int> &jump_pair);

  private:
    // Parameters info
    Config config_;
    long long config_dump_{};
    long long log_dump_{};
    double first_cutoff_ = Al_const::kFirstNearestNeighborCutoff;
    double second_cutoff_ = Al_const::kSecondNearestNeighborsCutoff;
    // Statistical info
    double time_{};
    double total_energy_{};
    long long iter_{};
    long long step_{};

    mutable boost::mpi::environment env_;
    mutable boost::mpi::communicator world_;
    int mpi_rank_;
    int mpi_size_;
    mutable std::mt19937_64 generator_;
    mutable std::ofstream ofs_;
};

} // namespace kn
#endif //KN_INCLUDE_KMCSIMULATION_H_
