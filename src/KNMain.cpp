#include "ConfigModifier.h"
using namespace box;
using namespace std;
int main(int argc, char *argv[]) {
  Config test;
  // test.GenerateFCC(4.046,"Al",{30,30,30});
  // test.GenerateHCP(3.209,5.211,"Mg",{50,50,50});

  // if (!test.ReadPOSCAR("test.pos")) { return 5; }

  // test.Perturb();
  test.ReadConfig("0.cfg");


  // test.ConvertRelativeToAbsolute();
  // test.ConvertAbsoluteToRelative();

  // test.UpdateNeighbors(3.5,4.5);
  test.ShiftAtomToCentral(0);
  Vector3<double> plane{1, 1, 1};
  double plane_inner_product =  InnerProduct(plane);
  Vector3<double> planeMoveDistance = plane*(4.046/plane_inner_product/2);
  test.MoveAbsoluteDistance(planeMoveDistance);

  map<std::string, int> bond = test.CountAllBonds(3.5);

  for (const auto &[element, count]:bond) {
    std::cout << element<<": ";
    std::cout << count << "\n";
  }

  test.WriteConfig();
  test.WritePOSCAR();

  return 0;
}
