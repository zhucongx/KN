#include "MpiIterator.h"
MpiIterator::MpiIterator() {}
MpiIterator::MpiIterator(long long int initial_number,
                         long long int increment_number,
                         long long int finial_number)
    : initial_number_(initial_number),
      increment_number_(increment_number),
      finial_number_(finial_number) {

  mpi_rank_ = world_.rank();
  mpi_size_ = world_.size();
}
void MpiIterator::IterateToFindCLusters() {
  long long num_total_loop = (finial_number_ - initial_number_) / increment_number_ + 1;
  long long quotient = num_total_loop / mpi_size_;
  int remainder = num_total_loop % mpi_size_;
  long long num_cycle = remainder ? (quotient + 1) : quotient;
  for (long long i = 0; i < num_cycle; i++) {
    long long num_file = initial_number_ + (i * quotient + mpi_rank_) * increment_number_;
    if (num_file > num_total_loop)
      break;
    std::string file_name = std::to_string(num_file) + ".cfg";
    box::ClusterFinder cluster_finder(file_name, "Al", 3, 3);
    cluster_finder.FindClustersAndOutput();
  }
  
}
