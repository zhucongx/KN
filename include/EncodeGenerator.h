#ifndef KN_INCLUDE_ENCODEGENERATOR_H_
#define KN_INCLUDE_ENCODEGENERATOR_H_

#include <vector>
#include <map>
#include "Bond.h"
#include "Config.h"

namespace kn {
class EncodeGenerator {
  public:
    EncodeGenerator(std::string reference_filename);
    void PrintOutEncode();
    static std::vector<std::array<std::string, Al_const::kLengthOfEncodes>> Encode(
        const Config &config,
        const std::pair<int, int> &jump_pair);
  private:
    std::string reference_filename_;
};
} // namespace kn
#endif //KN_INCLUDE_ENCODEGENERATOR_H_
