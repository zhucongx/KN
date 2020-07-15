#ifndef KN_INCLUDE_ENCODEGENERATOR_H_
#define KN_INCLUDE_ENCODEGENERATOR_H_

#include <vector>
#include <map>
#include "Bond.h"
#include "Config.h"

namespace kn ::EncodeGenerator {
// 0. Atom2
//
// 1. Atom1.First && Atom2.First 4
// 2. Atom1.First && Atom2.Second || Atom1.Second && Atom2.First 2+2
// 3. Atom1.First && Atom2.Third || Atom1.Third && Atom2.First 4+4
// 4. Atom1.First && not Atom2.NNL || Atom2.First && not Atom1.NLL 2
//
// 5. Atom1.Second && Atom2.Second 0
// 6. Atom1.Second && Atom2.Third || Atom1.Third && Atom2.Second 2+2
// 7. Atom1.Second && not Atom2.NNL || Atom2.Second && not Atom1.NLL 2+2
//
// 8. Atom1.Third && Atom2.Third 4
// 9. Atom1.Third && not Atom2.NNL || Atom2.Third && not Atom1.NLL 14+14

constexpr int kNumOfSubEncode = 10;
enum LengthOfEncodes {
  kLengthOfZerothEncodes = 1,
  kLengthOfFirstEncodes = 4,
  kLengthOfSecondEncodes = 4,
  kLengthOfThirdEncodes = 8,
  kLengthOfFourthEncodes = 2,
  kLengthOfFifthEncodes = 0,
  kLengthOfSixthEncodes = 4,
  kLengthOfSeventhEncodes = 4,
  kLengthOfEighthEncodes = 4,
  kLengthOfNinthEncodes = 28,
};
const std::array<LengthOfEncodes, kNumOfSubEncode> All_Encode_Length{kLengthOfZerothEncodes,
                                                                     kLengthOfFirstEncodes,
                                                                     kLengthOfSecondEncodes,
                                                                     kLengthOfThirdEncodes,
                                                                     kLengthOfFourthEncodes,
                                                                     kLengthOfFifthEncodes,
                                                                     kLengthOfSixthEncodes,
                                                                     kLengthOfSeventhEncodes,
                                                                     kLengthOfEighthEncodes,
                                                                     kLengthOfNinthEncodes};

constexpr int kLengthOfEncodes =
    kLengthOfZerothEncodes + kLengthOfFirstEncodes + kLengthOfSecondEncodes + kLengthOfThirdEncodes
        + kLengthOfFourthEncodes + kLengthOfFifthEncodes + kLengthOfSixthEncodes
        + kLengthOfSeventhEncodes + kLengthOfEighthEncodes + kLengthOfNinthEncodes;

std::vector<std::array<int, kLengthOfEncodes>> Encode(
    const Config &config,
    const std::pair<int, int> &jump_pair,
    const std::unordered_map<std::string, int> &type_category_hashmap);

std::array<int, kLengthOfEncodes> GetBackwardEncode(
    const std::array<int, kLengthOfEncodes> &forward_encode);

void PrintOutEncode(const std::string &reference_filename,
                    const std::unordered_map<std::string, int> &type_category_hashmap);
} // namespace kn
#endif //KN_INCLUDE_ENCODEGENERATOR_H_
