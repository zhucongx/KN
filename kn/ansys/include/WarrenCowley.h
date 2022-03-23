//
// Created by Zhucong Xi on 3/15/22.
//

#ifndef KN_KN_ANSYS_INCLUDE_WARRENCOWLEY_H_
#define KN_KN_ANSYS_INCLUDE_WARRENCOWLEY_H_

#include <vector>
#include <set>
#include <unordered_set>

#include "Config.h"

namespace ansys {

class WarrenCowley {
  public:
    explicit WarrenCowley(const std::string &cfg_filename);
    [[nodiscard]] std::map<std::string, double> FindWarrenCowley() const;
  protected:
    cfg::Config config_;
    std::set<std::string> element_set_{};

};

} // namespace ansys


#endif //KN_KN_ANSYS_INCLUDE_WARRENCOWLEY_H_
