#ifndef KN_KN_KMC_INCLUDE_DFTANALYSIS_H_
#define KN_KN_KMC_INCLUDE_DFTANALYSIS_H_

#include "string"
#include <unordered_map>
#include <unordered_set>
// This namespace contains some functions that can analyze dft calculations data
namespace kn::DftAnalysis {

// Generate encode file according to the log.txt file and generated config files
void PrintOutEncode(const std::string &reference_filename,
                    const std::unordered_map<std::string, int> &type_category_hashmap);

void PrintOutClusterExpansionAverage(
    const std::string &reference_filename,
    const std::unordered_map<std::string, double> &type_category_hashmap);


// // Generate bond change file according to the log.txt file and generated config files
// void PrintOutBond(const std::string &reference_name,
//                   const std::unordered_set<std::string> &type_hashset);
} // namespace kn::DftAnalysis

#endif //KN_KN_KMC_INCLUDE_DFTANALYSIS_H_
