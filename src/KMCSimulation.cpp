#include "KMCSimulation.h"
#include <chrono>

#include "ConfigIO.h"
#include "EncodeGenerator.h"
namespace kn {

KMCSimulation::KMCSimulation(const std::string &cfg_filename,
                             double first_nearest_neighbors_distance,
                             double second_nearest_neighbors_distance) :
    generator_(std::chrono::system_clock::now().time_since_epoch().count()) {
  config_ = ConfigIO::ReadConfig(cfg_filename, true);
  mpi_rank_ = world_.rank();
  mpi_size_ = world_.size();
  if (mpi_rank_ == 0) {
    std::cout << "Using " << mpi_size_ << " processes." << std::endl;
    ofs_.open("kmc.txt", std::ofstream::out | std::ofstream::app);
  }
}
void KMCSimulation::RunSimulation() {
  std::uniform_real_distribution<> dis(0, 1);
}

std::vector<double> KMCSimulation::CalculateBarrierAndEnergyDifference(
    const std::pair<int, int> &jump_pair) const {
  const auto[first, second] = jump_pair;
  std::vector<std::string> codes; // atom location in the original atom list
  const auto input_forward = EncodeGenerator::Encode(config_, jump_pair);
  int num_row = input_forward.size();

  // const auto &  = encodes;


}

} // namespace kn
