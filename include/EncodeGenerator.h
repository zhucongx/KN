#ifndef KN_INCLUDE_ENCODEGENERATOR_H_
#define KN_INCLUDE_ENCODEGENERATOR_H_

#include <vector>
#include <map>
#include "Bond.h"
#include "Config.h"

namespace kn :: EncodeGenerator {
// public:
void PrintOutEncode(const std::string &reference_filename);
std::vector<std::array<std::string, Al_const::kLengthOfEncodes>> Encode(
    const Config &config,
    const std::pair<int, int> &jump_pair);

};
} // namespace kn
#endif //KN_INCLUDE_ENCODEGENERATOR_H_
