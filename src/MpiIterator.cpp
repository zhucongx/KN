#include "MpiIterator.h"
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
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
void MpiIterator::IterateToFindCLusters(const std::string &solvent_atom_type,
                                        int smallest_cluster_criteria,
                                        int solvent_bond_criteria) {
  long long num_total_loop = (finial_number_ - initial_number_) / increment_number_ + 1;
  long long quotient = num_total_loop / mpi_size_;
  auto remainder = num_total_loop % mpi_size_;
  long long num_cycle = remainder ? (quotient + 1) : quotient;

  for (long long i = 0; i < num_cycle; i++) {
    long long num_file = initial_number_ + (i * mpi_size_ + mpi_rank_) * increment_number_;
    if (num_file > finial_number_)
      break;
    std::string file_name = std::to_string(num_file) + ".cfg";

    box::ClusterFinder cluster_finder(file_name,
                                      solvent_atom_type,
                                      smallest_cluster_criteria,
                                      solvent_bond_criteria);
    auto num_different_element = cluster_finder.FindClustersAndOutput();
    world_.barrier();
    if (mpi_rank_ != 0) {
      boost::mpi::gather(world_, num_different_element, 0);
    } else {
      std::vector<box::ClusterFinder::ClusterElementNumMap> all_num_different_element;
      boost::mpi::gather(world_, num_different_element, all_num_different_element, 0);
      std::ofstream ofs("clusters_info.txt", std::ofstream::out | std::ofstream::app);
      auto file_index = num_file;
      for (const auto &this_num_different_element:all_num_different_element) {
        ofs << "Config " << file_index << '\n';
        for (const auto &cluster:this_num_different_element) {
          for(const auto&[key,count]:cluster){
            ofs << count << ' ';
          }
          ofs << '\n';
        }
        file_index += increment_number_;
      }

    }
  }
}
