#ifndef KN_INCLUDE_CONFIGIO_H_
#define KN_INCLUDE_CONFIGIO_H_
#include "Config.h"
namespace kn {
class ConfigIO {
 public:
  static Config ReadPOSCAR(const std::string &filename);
  static Config ReadConfig(const std::string &filename);

  // Write Configuration out as POSCAR file. If the show_vacancy_option is
  // true, output will have "X" for visualization. If false, vacancies will be
  // ignored for VASP calculation.
  static void WritePOSCAR(const Config &config,
                          const std::string &filename,
                          bool show_vacancy_option = false);
  static void WriteConfig(const Config &config,
                          const std::string &filename,
                          bool neighbors_info = true);
};

}// namespace kn
#endif //KN_INCLUDE_CONFIGIO_H_
