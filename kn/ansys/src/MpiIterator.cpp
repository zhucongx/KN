#include "MpiIterator.h"
#include <iostream>
#include <mpi.h>

namespace ansys {

MpiIterator::MpiIterator(unsigned long long int initial_number,
                         unsigned long long int increment_number,
                         unsigned long long int finial_number)
    : initial_number_(initial_number),
      increment_number_(increment_number),
      finial_number_(finial_number) {
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
  if (mpi_rank_ == 0) {
    std::cout << "Using " << mpi_size_ << " processes." <<
              std::endl;
  }
}
MpiIterator::~MpiIterator() {
  MPI_Finalize();
}
} // namespace ansys
