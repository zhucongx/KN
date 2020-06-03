#include "MpiClusters.h"
#include <utility>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

namespace kn {
MpiClusters::MpiClusters(long long int initial_number,
                         long long int increment_number,
                         long long int finial_number,
                         std::string solvent_element,
                         int smallest_cluster_criteria,
                         int solvent_bond_criteria) :
  MpiIterator(initial_number,
              increment_number,
              finial_number),
  solvent_element_(std::move(solvent_element)),
  smallest_cluster_criteria_(smallest_cluster_criteria),
  solvent_bond_criteria_(solvent_bond_criteria) {
}

void MpiClusters::IterateToRun() {
  long long num_total_loop = (finial_number_ - initial_number_) / increment_number_ + 1;
  long long quotient = num_total_loop / mpi_size_;
  auto remainder = num_total_loop % mpi_size_;
  long long num_cycle = remainder ? (quotient + 1) : quotient;

  for (long long i = 0; i < num_cycle; ++i) {
    long long num_file = initial_number_ + (i * mpi_size_ + mpi_rank_) * increment_number_;
    ClustersFinder::ClusterElementNumMap num_different_element;
    if (num_file <= finial_number_) {
      std::string filename = std::to_string(num_file) + ".cfg";

      ClustersFinder cluster_finder(filename,
                                    solvent_element_,
                                    smallest_cluster_criteria_,
                                    solvent_bond_criteria_);
      num_different_element = cluster_finder.FindClustersAndOutput();
    }
    world_.barrier();
    if (mpi_rank_ != 0) {
      boost::mpi::gather(world_, num_different_element, 0);
    } else {
      std::vector<ClustersFinder::ClusterElementNumMap> all_num_different_element;
      boost::mpi::gather(world_, num_different_element, all_num_different_element, 0);

      std::ofstream ofs("clusters_info.txt", std::ofstream::out | std::ofstream::app);
      auto file_index = num_file;
      for (const auto &this_num_different_element : all_num_different_element) {
        if (file_index > finial_number_)
          break;
        ClustersFinder::PrintLog(std::to_string(file_index), this_num_different_element);
        file_index += increment_number_;
      }
    }
  }
}
} // namespace kn
