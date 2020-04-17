#include "BondCounter.h"
using namespace box;
using namespace std;

int main(int argc, char *argv[]) {
  long long start = 1579000000;
  long long interval = 1000000;
  long long end = 2256000000;
  map<Bond, double> bond_energy{{{"Al", "Al"}, {-0.60648589}},
                                {{"Al", "Mg"}, {-0.41880444}},
                                {{"Al", "Zn"}, {-0.38473348}},
                                {{"Mg", "Mg"}, {-0.21746198}},
                                {{"Mg", "Zn"}, {-0.23373609}},
                                {{"Zn", "Zn"}, {-0.17651597}}};
  // start = 1;
  // end = 1;
  for (long long i = start; i <= end; i += interval) {
    double energy = 0;
    string fname = to_string(i);
    fname += ".cfg";
    BondCounter test({30, 30, 30}, {1, 1, 1}, {0.5, 0.5, 0});
    test.config_.ReadConfig(fname);
    std::map<Bond, int> bonds_changed = test.GetBondChange();

    cout << "#" << fname << endl;
    for (const auto&[key, count] : bonds_changed) {
      cout << "#" << key << " " << count << '\n';
      energy += static_cast<double>(count) * bond_energy[key];
    }
    cout << energy << '\n';
  }

  // Config test;
  // test.GenerateFCC(4.046,"Al",{2,2,2});
  // test.MoveRelativeDistance({1./6.,1/6.,1./6.});
  // test.WriteConfig();


  return 0;
}
