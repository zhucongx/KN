#ifndef KN_KN_KMC_INCLUDE_MPIITERATOR_H_
#define KN_KN_KMC_INCLUDE_MPIITERATOR_H_

#include <boost/mpi.hpp>
#include <ClustersFinder.h>

namespace kn {

class MpiIterator {
  public:
    virtual void IterateToRun() = 0;

    MpiIterator(long long int initial_number,
                long long int increment_number,
                long long int finial_number);
  protected:
    boost::mpi::environment env_;
    boost::mpi::communicator world_;
    int mpi_rank_;
    int mpi_size_;

    long long initial_number_;
    long long increment_number_;
    long long finial_number_;
};

} // namespace kn


#endif //KN_KN_KMC_INCLUDE_MPIITERATOR_H_
