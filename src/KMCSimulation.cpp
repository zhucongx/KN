#include "KMCSimulation.h"
#include <chrono>

#include <tensorflow/c/c_api.h>
#include "Model.h"
#include "Tensor.h"

#include "EncodeGenerator.h"
#include "LRUCache.h"
namespace kn {

KMCSimulation::KMCSimulation(const std::string &cfg_filename) :
    generator_(std::chrono::system_clock::now().time_since_epoch().count()) {
  config_ = Config::ReadConfig(cfg_filename, true);
  if (world_.rank() == 0) {
    std::cout << "Using " << world_.size() << " processes." << std::endl;
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
