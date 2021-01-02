#ifndef KN_KN_API_INCLUDE_PARAMETER_H_
#define KN_KN_API_INCLUDE_PARAMETER_H_
#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace api {
class Parameter {
  public:
    void ParseArgs(int argc, char* argv[]) {
      for (int i = 0; i < argc; ++i) {
        if (!strcmp((argv[i]), "--p") || !strcmp(argv[i], "-p"))
          sparams["parameter_file"] = string(argv[++i]);
        if (!strcmp(argv[i], "--i") || !strcmp(argv[i], "-i"))
          sparams["dmpfile"] = string(argv[++i]);
        if (!strcmp(argv[i], "--f") || !strcmp(argv[i], "-f"))
          sparams["potfile"] = string(argv[++i]);
      }
    }

};


}

#endif //KN_KN_API_INCLUDE_PARAMETER_H_
