#ifndef KN_INCLUDE_CONFIGMODIFIER_H_
#define KN_INCLUDE_CONFIGMODIFIER_H_

#include "Cfg.hpp"
#include "FCCConfig.h"

class ConfigModifier {
  public:

    static ConfigModifier *Instance();

    void RegisterConfig(int index, box::Config *registrar);
    box::Config *GetConfig(int index);
    int GetNumConfig();
    // forbid creating from outside
    ConfigModifier() = delete;
    ~ConfigModifier() = delete;
    // forbid copying and assigning
    ConfigModifier(const ConfigModifier &) = delete;
    const ConfigModifier &operator=(const ConfigModifier &) = delete;
  private:

    static ConfigModifier *instance;
    std::map<int, box::Config*> map_config_;

};

ConfigModifier *ConfigModifier::instance = new(std::nothrow) ConfigModifier;

ConfigModifier *ConfigModifier::Instance() {
  return instance;
}

void ConfigModifier::RegisterConfig(int index, box::Config *registrar) {
  map_config_[index] = registrar;
}

box::Config *ConfigModifier::GetConfig(int index) {
  auto it = map_config_.find(index);
  if (it != map_config_.end()) {
    return it->second;
  }
  std::cout << "error";
  return nullptr;
}

int ConfigModifier::GetNumConfig() {
  return map_config_.size();
}

ConfigModifier::ConfigModifier() {
}

#endif //KN_INCLUDE_CONFIGMODIFIER_H_
