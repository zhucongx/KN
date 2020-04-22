#include "BondCounter.h"
using namespace box;
using namespace std;

int main(int argc, char *argv[]) {
  long long start = 0;
  long long interval = 1000000;
  long long end = 264000000;
  // end = 0;
  // map<Bond, double> bond_energy{{{"Al", "Al"}, {-0.60648589}},
  //                               {{"Al", "Mg"}, {-0.41880444}},
  //                               {{"Al", "Zn"}, {-0.38473348}},
  //                               {{"Mg", "Mg"}, {-0.21746198}},
  //                               {{"Mg", "Zn"}, {-0.23373609}},
  //                               {{"Zn", "Zn"}, {-0.17651597}}};
  map<Bond, map<long long, int>> bond_store;
  std::ofstream ofs("COUTPUT.txt", std::ofstream::out);
  for (long long i = start; i <= end; i += interval) {
    string fname = to_string(i);
    fname += ".cfg";
    BondCounter test({30, 30, 30}, {1, 1, 1}, {0.5, 0.5, 0});
    // BondCounter test;
    // test.factor_ = {30, 30, 30};
    // test.plane_set_.insert({-1, -1, -1});
    // test.SetBurgersVector({0.5, 0.5, 0});
    test.config_.ReadConfig(fname);
    std::map<Bond, int> bonds_changed = test.GetBondChange();

    ofs << "#" << fname << endl;
    for (const auto&[key, count] : bonds_changed) {
      ofs << "#" << key << " " << count << '\n';
      bond_store[key][i] = count;
      // energy += static_cast<double>(count) * bond_energy[key];
    }

    ofs << '\n';

    for (const auto&[key, map_count] : bond_store) {
      ofs << key << ": ";
      for (long long j = start; j <= i; j += interval) {
        auto it = map_count.find(j);
        ofs << (it == map_count.end() ? 0 : it->second) << " ";
      }
      ofs << '\n';
      // energy += static_cast<double>(count) * bond_energy[key];
    }

    ofs << '\n';
  }

  // Config test;
  // test.GenerateHCP(4, 9, "Al", {30, 30, 30});
  // test.ConvertAbsoluteToRelative();
  // test.ConvertRelativeToAbsolute();
  //
  // // test.MoveRelativeDistance({1./6.,1/6.,1./6.});
  // test.WriteConfig();


  return 0;
}
