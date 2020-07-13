#include "MpiClusters.h"
#include "ConfigIO.h"
#include "EncodeGenerator.h"
#include "KMCSimulation.h"
#include "ConfigGenerator.h"
using namespace kn;
using namespace std;
namespace mpi = boost::mpi;

int main(int argc, char *argv[]) {

  // EncodeGenerator dw("log.txt");
  // dw.PrintOutEncode();

  Config test = ConfigIO::ReadConfig("0.cfg", true);

  auto asd = EncodeGenerator::Encode(test, {82, 83});

  for (auto &&a:asd) {
    for (auto &&b:a)
      cout << b << "   ";
    cout << '\n';
  }

  auto asd1 = EncodeGenerator::Encode(test, {83, 82});
  for (auto &&a:asd1) {
    for (auto it = a.rbegin(); it != a.rend(); it ++)
      cout << *it << "   ";
    cout << '\n';
  }


  // test.UpdateNeighbors(Al_const::kFirstNearestNeighborsCutoff,
  //                      Al_const::kNearNeighborsCutoff);
  // ConfigIO::WriteConfig(test, "0_f.cfg", false);
  // MpiClusters test(370000000, 1000000, 6252000000,
  //                  "Al", 3, 3);
  // test.IterateToRun();

  // MpiNeighbors test(370000000, 1000000, 6252000000);
  // test.IterateToRun();

  // long long start = 0;
  // long long interval = 1000000;
  // long long end = 264000000;
  // // end = 0;
  // // map<Bond, double> bond_energy{{{"Al", "Al"}, {-0.60648589}},
  // //                               {{"Al", "Mg"}, {-0.41880444}},
  // //                               {{"Al", "Zn"}, {-0.38473348}},
  // //                               {{"Mg", "Mg"}, {-0.21746198}},
  // //                               {{"Mg", "Zn"}, {-0.23373609}},
  // //                               {{"Zn", "Zn"}, {-0.17651597}}};
  // map<Bond, map<long long, int>> bond_store;
  // std::ofstream ofs("COUTPUT.txt", std::ofstream::out);
  // for (long long i = start; i <= end; i += interval)
  // {
  //   string fname = to_string(i);
  //   fname += ".cfg";
  //   BondCounter test(fname,{30, 30, 30}, {1, 1, 1}, {0.5, 0.5, 0});
  //   // BondCounter test;
  //   // test.factor_ = {30, 30, 30};
  //   // test.plane_set_.insert({-1, -1, -1});
  //   // test.SetBurgersVector({0.5, 0.5, 0});
  //   std::map<Bond, int> bonds_changed = test.GetBondChange();
  //
  //   ofs << "#" << fname << endl;
  //   for (const auto&[key, count] : bonds_changed)
  //   {
  //     ofs << "#" << key << ' ' << count << '\n';
  //     bond_store[key][i] = count;
  //   }
  //
  //   ofs << '\n';
  //
  //   for (const auto&[key, map_count] : bond_store)
  //   {
  //     ofs << key << ": ";
  //     for (long long j = start; j <= i; j += interval)
  //     {
  //       auto it = map_count.find(j);
  //       ofs << (it == map_count.end() ? 0 : it->second) << ' ';
  //     }
  //     ofs << '\n';
  //     // energy += static_cast<double>(count) * bond_energy[key];
  //   }
  //
  //   ofs << '\n';
  // }

  //  BondCounter test({2, 2, 2}, {1, 1, 1}, {0.5, 0.5, 0});
  //  AntiPhaseConfig L10;
  //  L10.GenerateL10(4.046, {"Al", "Zn"}, {2, 2, 2});
  //  test.config_ = L10;
  //  L10.WritePOSCAR("0.poscar");
  //
  //  std::map<Bond, int> bonds_changed = test.GetBondChange();
  //
  //  std::cout << "#" << 0 << endl;
  //  for (const auto&[key, count] : bonds_changed)
  //  {
  //    std::cout << "#" << key << ' ' << count << '\n';
  //  }

  // Config cfg;
  // cfg.ReadPOSCAR("0.pos");
  // cfg.WritePOSCAR("00.pos");

  // Config cfg;
  // cfg.GenerateBCC(3.1652, "W", {10, 10, 10});
  // // cfg.LinearTransform({{0.57735,0.57735,0.57735},{0.408248,0.408248,-0.816497},{0.707107,-0.707107,0}});
  // cfg.WriteConfig("BCC.cfg");
  // return 0;
}
