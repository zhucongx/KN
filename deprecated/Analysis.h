#ifndef KN_KN_ANSYS_INCLUDE_ANALYSIS_H_
#define KN_KN_ANSYS_INCLUDE_ANALYSIS_H_
#include "string"
#include <unordered_map>

namespace ansys {
void PrintOutClusterExpansionAverage(
    const std::string &reference_filename,
    const std::unordered_map<std::string, double> &type_category_hashmap);


} // namespace ansys
#endif //KN_KN_ANSYS_INCLUDE_ANALYSIS_H_
