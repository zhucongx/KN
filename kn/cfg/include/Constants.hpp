#ifndef KN_KN_CFG_INCLUDE_CONSTANTS_HPP_
#define KN_KN_CFG_INCLUDE_CONSTANTS_HPP_

#include <string>
#include <vector>
#include <algorithm>

namespace Al_const {
constexpr double kFirstNearestNeighborsCutoff = 3.5;
constexpr double kSecondNearestNeighborsCutoff = 4.5;
constexpr double kThirdNearestNeighborsCutoff = 5.3;
// length between fourth nearest neighbors is double of length of first nearest neighbors
constexpr double kFourthNearestNeighborsCutoff = 5.9;
// constexpr double kNearNeighborsCutoff = kThirdNearestNeighborsCutoff;



// For FCC the first nearest neighbor is 12
constexpr size_t kNumFirstNearestNeighbors = 12;
// For FCC the second nearest neighbor is 6
constexpr size_t kNumSecondNearestNeighbors = 6;
// For FCC the third nearest neighbor is 6
constexpr size_t kNumThirdNearestNeighbors = 24;
// For FCC the fourth nearest neighbor is 12
constexpr size_t kNumFourthNearestNeighbors = 12;
} // namespace Al_const
#endif //KN_KN_CFG_INCLUDE_CONSTANTS_HPP_

