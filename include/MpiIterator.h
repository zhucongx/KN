#ifndef KN_INCLUDE_MPIIterator_H_
#define KN_INCLUDE_MPIIterator_H_

#include <boost/mpi.hpp>
#include <ClusterFinder.h>

class MpiIterator {
 public:
  MpiIterator(long long int initial_number,
              long long int increment_number,
              long long int finial_number);
  MpiIterator();

  void IterateToFindCLusters();
 private:
  boost::mpi::environment env_;
  boost::mpi::communicator world_;
  int mpi_rank_, mpi_size_;
  long long initial_number_;
  long long increment_number_;
  long long finial_number_;

};

#endif //KN_INCLUDE_MPIIterator_H_
