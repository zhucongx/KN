#include "Analysis.h"
#include <utility>
#include <mpi.h>
// #include <boost/serialization/string.hpp>
// #include <boost/serialization/vector.hpp>
// #include <boost/serialization/map.hpp>

#include "ClustersFinder.h"
#include "WarrenCowley.h"

namespace ansys {
Analysis::Analysis(unsigned long long int initial_number,
                   unsigned long long int increment_number,
                   std::string solvent_element,
                   size_t smallest_cluster_criteria,
                   size_t solvent_bond_criteria) :
    initial_number_(initial_number),
    increment_number_(increment_number),
    finial_number_(increment_number),
    solvent_element_(std::move(solvent_element)),
    smallest_cluster_criteria_(smallest_cluster_criteria),
    solvent_bond_criteria_(solvent_bond_criteria) {
  std::ifstream ifs("kmc_log_modified.txt", std::ifstream::in);
  if (!ifs.is_open()) {
    std::cout << "Cannot open kmc_log_modified.txt\n";
    return;
  }
  unsigned long long filename;
  double time;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  while (true) {
    if (ifs.eof() || ifs.bad()) {
      break;
    }
    if (ifs.peek() == '#') {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }

    ifs >> filename >> time;
    ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if (filename >= initial_number_ && (filename - initial_number_) % increment_number == 0) {
      filename_time_hashset_[filename] = time;
      finial_number_ = filename;
    }
  }
  finial_number_ -= increment_number;
}

Analysis::~Analysis() = default;

void Analysis::SerialRunCluster() const {
  // start
  std::ofstream ofs("clusters_info.json", std::ofstream::out);
  ofs << "[ \n";

  for (unsigned long long i = 0; i <= finial_number_; i += increment_number_) {
    ClustersFinder cluster_finder(std::to_string(i) + ".cfg",
                                  solvent_element_,
                                  smallest_cluster_criteria_,
                                  solvent_bond_criteria_);
    auto num_different_element = cluster_finder.FindClustersAndOutput();

    ofs << "{ \n"
        << "\"index\" : " << "\"" << std::to_string(i) << "\",\n"
        << "\"time\" : " << filename_time_hashset_.at(i) << ",\n"
        << "\"clusters\" : [ \n";
    for (auto it = num_different_element.cbegin(); it < num_different_element.cend(); ++it) {
      ofs << "[ ";
      const auto &cluster = *it;
      std::for_each(it->cbegin(), --(cluster.cend()), [&ofs](auto it) {
        ofs << it.second << ", ";
      });
      ofs << (cluster.crbegin())->second;
      if (it == num_different_element.cend() - 1) {
        ofs << "] \n";
      } else {
        ofs << "], \n";
      }
    }
    ofs << "]\n";
    if (i != finial_number_) {
      ofs << "}, \n";
    } else {
      ofs << "} \n";
    }
  }
  ofs << " ]" << std::endl;;
}

void Analysis::SerialRunWarrenCowley() const {
  // start
  std::ofstream ofs("warren_cowley.csv", std::ofstream::out);

  for (unsigned long long i = 0; i <= finial_number_; i += increment_number_) {
    WarrenCowley warren_cowley_finder(std::to_string(i) + ".cfg");
    if (i == 0) {
      for (const auto &[pair, value]: warren_cowley_finder.FindWarrenCowley()) {
        ofs << pair << "\t";
      }
      ofs << '\n';
    }
    ofs << i << "\t";
    for (const auto &[pair, value]: warren_cowley_finder.FindWarrenCowley()) {
      ofs << value << "\t";
    }
    ofs << std::endl;
  }
}
} // namespace ansys
