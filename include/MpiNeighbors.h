#ifndef KN_INCLUDE_MPINEIGHBORS_H_
#define KN_INCLUDE_MPINEIGHBORS_H_
#include "MpiIterator.h"
#include "ConfigIO.h"
namespace kn {

class MpiNeighbors : public MpiIterator {
  public:
    MpiNeighbors(long long int initial_number,
                 long long int increment_number,
                 long long int finial_number);
    void IterateToRun() override;
};

} // namespace kn
#endif //KN_INCLUDE_MPINEIGHBORS_H_
