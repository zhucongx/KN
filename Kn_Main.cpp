// #include "Encode.h"
// #include "ConfigGenerator.h"
#include <SizeMisfitGenerator.h>
#include "deprecated/Analysis.h"
// #include "ClusterExpansion.h"
// #include "NebDistance.h"
// using namespace kn;
using namespace std;
// namespace mpi = boost::mpi;
// pure Al lattice constant 4.0404778433873281
#include "KMCSimulation.h"
#include "MpiClusters.h"
// #include "ClusterExpansion.h"
int main(int argc, char *argv[]) {
  // gen::SizeMisfitGenerator a(4.0404778433873281,
  //                            {2, 2, 2},
  //                            "Al",
  //                            {"Al", "Cu", "Mg", "Zn", "X"},
  //                            "/Users/zhucongx/Program/goali/pot_old/potpaw_PBE/elements/");
  // a.CreateConfigs();
  // cfg::Config::WriteConfig(cfg::Config::ReadConfig("init.cfg"), "start.cfg");
  // ansys::MpiClusters test(0, 10000000, 180000000,
  //                         "Al", 3, 3);
  // test.IterateToRun();
  kmc::KMCSimulation a(cfg::Config::ReadConfig("start.cfg"),
                       1e5,
                       1e7,
                       1e10,
                       {"Al", "Mg", "Zn"},
                       0, 0, 0, "kmc_parameters.json", 3000);
  a.Simulate();

  // std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());
  // auto config = cfg::GenerateFCC(4.046, "Al", {30, 30, 30});
  // vector<size_t> index_v;
  // index_v.resize(108000);
  // for (size_t i = 0; i < 108000; ++i)
  //   index_v[i] = i;
  // shuffle(index_v.begin(), index_v.end(), generator);
  // for (size_t i = 0; i < 3240; ++i){
  //   config.atom_list_[index_v[i]].SetType("Mg");
  // }
  // for (size_t i = 3240; i < 3240+2700; ++i){
  //   config.atom_list_[index_v[i]].SetType("Zn");
  // }
  // cfg::Config::WriteConfig(config, "start_mf.cfg");
  // unordered_map<std::string, double> type_category_hashmap{{"Al", 0},
  //                                                          {"Mg", 2},
  //                                                          {"Zn", -1},
  //                                                          {"X",  3}};
  // ansys::PrintOutClusterExpansionAverage("log.txt", type_category_hashmap);

// auto cfg = ansys::ClusterExpansion::GetAverageClusterParametersForwardAndBackward(cfg::Config::ReadConfig("start.cfg"),
//                                                         {82, 83},
//                                                         type_category_hashmap)[0];
// for (auto kI : cfg) {
//   cout << kI << " ";
// }
// auto config = cfg::Config::ReadConfig("start.cfg");
// const auto move_distance = Vector_t{0.5, 0.5, 0.5} - GetPairCenter(config, {153,168});
// config.MoveRelativeDistance(move_distance);
// cfg::Config::WriteConfig(config,"new.cfg", false);


// NebDistance::PrintTheDistanceFromTwoPOSCARFiles("POSCAR0", "POSCAR1");


// test.UpdateNeighbors(Al_const::kFirstNearestNeighborsCutoff,
//                      Al_const::kNearNeighborsCutoff);
// ConfigIO::WriteConfig(test, "0_f.cfg", false);


// MpiNeighbors test(370000000, 1000000, 6252000000);
// test.IterateToRun();

// long long start = 0;

//     cout << b << " ";  // long long interval = 1000000;
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
// Config cfg;
// cfg.GenerateBCC(3.1652, "W", {10, 10, 10});
// // cfg.LinearTransform({{0.57735,0.57735,0.57735},{0.408248,0.408248,-0.816497},{0.707107,-0.707107,0}});
// cfg.WriteConfig("BCC.cfg");
// return 0;
}

//  std::cout << "#" << 0 << endl;
//  for (const auto&[key, count] : bonds_changed)
//  {
//    std::cout << "#" << key << ' ' << count << '\n';
//  }

// Config cfg;
// cfg.ReadPOSCAR("0.pos");
// cfg.WritePOSCAR("00.pos");
