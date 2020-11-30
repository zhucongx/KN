#include "MpiIterator.h"

namespace ansys {

MpiIterator::MpiIterator(unsigned long long int initial_number,
                         unsigned long long int increment_number,
                         unsigned long long int finial_number)
    : initial_number_(initial_number),
      increment_number_(increment_number),
      finial_number_(finial_number),
      mpi_rank_(world_.rank()),
      mpi_size_(world_.size()) {
  if (mpi_rank_ == 0) {
    std::cout << "Using " << mpi_size_ << " processes." <<
              std::endl;
  }
}
MpiIterator::~MpiIterator() = default;

} // namespace ansys