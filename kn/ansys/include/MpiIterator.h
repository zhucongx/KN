#ifndef KN_KN_KMC_INCLUDE_MPIITERATOR_H_
#define KN_KN_KMC_INCLUDE_MPIITERATOR_H_

namespace ansys {

class MpiIterator {
  public:
    virtual void IterateToRun() = 0;
    MpiIterator(unsigned long long int initial_number,
                unsigned long long int increment_number,
                unsigned long long int finial_number);
    virtual ~MpiIterator();

  protected:
    const unsigned long long initial_number_;
    const unsigned long long increment_number_;
    const unsigned long long finial_number_;

    int mpi_rank_;
    int mpi_size_;
};

} // namespace ansys


#endif //KN_KN_KMC_INCLUDE_MPIITERATOR_H_
