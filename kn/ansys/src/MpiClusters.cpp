#include "MpiClusters.h"
#include <utility>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include "ClustersFinder.h"

namespace ansys {
MpiClusters::MpiClusters(unsigned long long int initial_number,
                         unsigned long long int increment_number,
                         unsigned long long int finial_number,
                         std::string solvent_element,
                         size_t smallest_cluster_criteria,
                         size_t solvent_bond_criteria) :
    MpiIterator(initial_number,
                increment_number,
                finial_number),
    solvent_element_(std::move(solvent_element)),
    smallest_cluster_criteria_(smallest_cluster_criteria),
    solvent_bond_criteria_(solvent_bond_criteria) {
  if (mpi_rank_ == 0) {
    std::ifstream ifs("kmc_log.txt", std::ifstream::in);
    if (!ifs.is_open()) {
      std::cout << "Cannot open kmc_log.txt\n";
      return;
    }
    unsigned long long filename;
    double time;
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifs >> filename >> time) {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      if (filename >= initial_number_ && filename <= finial_number_
          && (filename - initial_number_) % increment_number == 0)
        filename_time_hashset_[filename] = time;
    }
    for (const auto &[element_type, index_vector] : filename_time_hashset_) {
      std::cout << element_type << " " << index_vector << "\n";
    }
  }
}

MpiClusters::~MpiClusters() = default;

void MpiClusters::IterateToRun() {
  auto num_total_loop = (finial_number_ - initial_number_) / increment_number_ + 1;
  auto quotient = num_total_loop / static_cast<unsigned long long int>(mpi_size_);
  auto remainder = num_total_loop % static_cast<unsigned long long int>(mpi_size_);
  auto num_cycle = remainder ? (quotient + 1) : quotient;
  // start
  std::ofstream ofs("clusters_info.json", std::ofstream::out);
  if (mpi_rank_ == 0) {
    ofs << "[ \n";
  }

  for (unsigned long long i = 0; i < num_cycle; ++i) {
    auto num_file = initial_number_ +
        (i * static_cast<unsigned long long int>(mpi_size_)
            + static_cast<unsigned long long int>(mpi_rank_)) * increment_number_;
    ClustersFinder::ClusterElementNumMap num_different_element;
    if (num_file <= finial_number_) {
      ClustersFinder cluster_finder(std::to_string(num_file) + ".cfg",
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

      auto file_index = num_file;
      for (const auto &this_num_different_element : all_num_different_element) {
        if (file_index > finial_number_)
          break;
        ofs << "{ \n"
            << "\"index\" : " << "\"" << std::to_string(file_index) << "\",\n"
            << "\"time\" : " << 10*filename_time_hashset_[file_index] << ",\n"
            << "\"clusters\" : [ \n";

        for (size_t i = 0; i < this_num_different_element.size(); i++) {
          const auto &cluster = this_num_different_element[i];
          ofs << "[ ";

          std::for_each(cluster.cbegin(), --cluster.cend(), [&ofs](auto it) {
            ofs << it.second << ", ";
          });
          ofs << (--cluster.cend())->second;
          if (i == this_num_different_element.size() - 1) {
            ofs << "] \n";
          } else {
            ofs << "], \n";
          }
        }
        ofs << "]\n";
        ofs << "}, \n";
        file_index += increment_number_;
      }
    }
  }

  if (mpi_rank_ == 0) {
    ofs << " ]";
  }
}
} // namespace ansys
