#ifndef KN_KN_KMC_INCLUDE_MPIITERATOR_H_
#define KN_KN_KMC_INCLUDE_MPIITERATOR_H_

#include <boost/mpi.hpp>
#include <ClustersFinder.h>

namespace kn {

class MpiIterator {
  public:
    virtual void IterateToRun() = 0;
    MpiIterator(unsigned long long int initial_number,
                unsigned long long int increment_number,
                unsigned long long int finial_number);
    virtual ~MpiIterator();

  protected:
    boost::mpi::environment env_;
    boost::mpi::communicator world_;
    int mpi_rank_;
    int mpi_size_;

    unsigned long long initial_number_;
    unsigned long long increment_number_;
    unsigned long long finial_number_;
};

} // namespace kn


#endif //KN_KN_KMC_INCLUDE_MPIITERATOR_H_
