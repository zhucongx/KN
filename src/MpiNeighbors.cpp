#include "MpiNeighbors.h"

namespace kn {

MpiNeighbors::MpiNeighbors(long long int initial_number,
                           long long int increment_number,
                           long long int finial_number)
    : MpiIterator(initial_number, increment_number, finial_number) {}
void MpiNeighbors::IterateToRun() {
  long long num_total_loop = (finial_number_ - initial_number_) / increment_number_ + 1;
  long long quotient = num_total_loop / mpi_size_;
  auto remainder = num_total_loop % mpi_size_;
  long long num_cycle = remainder ? (quotient + 1) : quotient;

  for (long long i = 0; i < num_cycle; ++i) {
    long long num_file = initial_number_ + (i * mpi_size_ + mpi_rank_) * increment_number_;
    ClustersFinder::ClusterElementNumMap num_different_element;
    if (num_file > finial_number_)
      break;

    auto config = ConfigIO::ReadConfig(std::to_string(num_file) + ".cfg", true);
    ConfigIO::WriteConfig(config, std::to_string(num_file) + "_neighbor.cfg", true);
  }
}
} // namespace kn
