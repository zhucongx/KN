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
    KMCSimulation(const std::string& cfg_filename);
    void RunSimulation();
    std::vector<double> CalculateBarrierAndEnergyDifference(const std::pair<int,
                                                                            int> &jump_pair) const;

  private:
    // Parameters info
    Config config_;
    long long config_dump_{};
    long long log_dump_{};

    // Statistical info
    double time_{};
    double total_energy_{};
    long long iter_{};
    long long step_{};

    // Running parameters
    mutable boost::mpi::environment env_;
    mutable boost::mpi::communicator world_;
    mutable std::mt19937_64 generator_;
    mutable std::ofstream ofs_;
};

} // namespace kn
#endif //KN_INCLUDE_KMCSIMULATION_H_
