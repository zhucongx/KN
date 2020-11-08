#include "Analysis.h"

#include <fstream>
#include <sstream>
#include <set>

#include "ClusterExpansion.h"
namespace ansys {
void PrintOutClusterExpansionAverage(
    const std::string &reference_filename,
    const std::unordered_map<std::string, double> &type_category_hashmap) {
  std::ifstream ifs(reference_filename, std::ifstream::in);
  std::ofstream ofs("cluster_expansion.txt", std::ofstream::out);
  std::string buffer;
  ofs << "config image ";
  for (size_t i = 0; i < 194; ++i) {
    ofs << "A" << i << " ";
  }
  for (size_t i = 0; i < 194; ++i) {
    ofs << "B" << i << " ";
  }
  ofs << '\n';

  while (ifs >> buffer) {
    if (buffer != "config") {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    int config_index, image_index, jump_pair_first, jump_pair_second;
    // config 0 end 0 pair: 248 218
    ifs >> config_index >> buffer >> image_index >> buffer >> jump_pair_first >> jump_pair_second;

    const auto
        config =
        cfg::Config::ReadConfig("config" + std::to_string(config_index) + "/s/start.cfg", true);
    auto cluster_expansion_average_code = ClusterExpansion::GetAverageClusterParameters(
        config, {jump_pair_first, jump_pair_second}, type_category_hashmap);
    auto cluster_expansion_average_code_back = ClusterExpansion::GetAverageClusterParametersBack(
        config, {jump_pair_first, jump_pair_second}, type_category_hashmap);

    ofs << config_index << "  " << image_index << "  ";
    for (const auto &code : cluster_expansion_average_code) {
      ofs << code << " ";
    }
    for (const auto &code : cluster_expansion_average_code_back) {
      ofs << code << " ";
    }
    ofs << '\n';
  }
}

} // namespace ansys