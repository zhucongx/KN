#include "DftAnalysis.h"

#include <fstream>
#include <sstream>
#include <set>
#include "Encode.h"
#include "ClusterExpansion.h"

namespace kn {
void DftAnalysis::PrintOutEncode(
    const std::string &reference_filename,
    const std::unordered_map<std::string, int> &type_category_hashmap) {
  std::ifstream ifs(reference_filename, std::ifstream::in);
  std::ofstream ofs("encode.txt", std::ofstream::out);
  std::string buffer;
  ofs << "config image ";
  for (int i = 0; i < Encode::kLengthOfEncodes; ++i) {
    ofs << "F" << i << " ";
  }
  for (int i = 0; i < Encode::kLengthOfEncodes; ++i) {
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

    auto config = Config::ReadConfig("config" + std::to_string(config_index) + "/s/start.cfg",
                                     true);
    auto forward_encode_result = Encode::GetEncode(
        config, {jump_pair_first, jump_pair_second}, type_category_hashmap);
    for (const auto &image_forward_encode : forward_encode_result) {
      ofs << config_index << "  " << image_index << "  ";
      for (const auto &code : image_forward_encode) {
        ofs << code << " ";
      }
      for (const auto &code : Encode::GetBackwardEncode(image_forward_encode)) {
        ofs << code << " ";
      }
      ofs << '\n';
    }
  }
}

// void DftAnalysis::PrintOutBond(const std::string &reference_name,
//                                const std::unordered_set<std::string> &type_hashset) {
//   std::ifstream ifs(reference_name, std::ifstream::in);
//   std::ofstream ofs("bond.txt", std::ofstream::out);
//   std::string buffer;
//   std::set<Bond> bonds_set;
//   for (const auto &element_type1 : type_hashset) {
//     for (const auto &element_type2 : type_hashset) {
//       bonds_set.insert({element_type1, element_type2});
//     }
//   }
//   ofs << "config image element ";
//   for (const auto &bond : bonds_set) {
//     ofs << bond << "_before ";
//   }
//   for (const auto &bond : bonds_set) {
//     ofs << bond << "_after ";
//   }
//   ofs << '\n';
//
//   while (ifs >> buffer) {
//     if (buffer != "config") {
//       ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//       continue;
//     }
//     int config_index, image_index, jump_pair_first, jump_pair_second;
//     // config 0 end 0 pair: 248 218
//     ifs >> config_index >> buffer >> image_index >> buffer >> jump_pair_first >> jump_pair_second;
//
//     auto config = Config::ReadConfig("config" + std::to_string(config_index) + "/s/start.cfg",
//                                      true);
//     auto[bond_map_before, bond_map_after] = Encode::GetBondAroundPair(
//         config, {jump_pair_first, jump_pair_second});
//
//     ofs << config_index << "  " << image_index << "  ";
//     ofs << config.GetAtomList()[jump_pair_second].GetType() << "  ";
//     for (const auto &bond : bonds_set) {
//       ofs << bond_map_before[bond] << " ";
//     }
//     for (const auto &bond : bonds_set) {
//       ofs << bond_map_after[bond] << " ";
//     }
//     ofs << '\n';
//   }
// }

} // namespace kn