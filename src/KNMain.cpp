#include "FCCConfig.h"
using namespace box;
int main(int argc, char *argv[]) {
  Config test;
  // test.GenerateFCC(4.046,"Al",{30,30,30});

  test.GenerateHCP(3.209,5.211,"Mg",{50,50,50});

  // if(!test.ReadPOSCAR("test.pos")) { return 5;}

  test.Perturb();
  // test.ReadConfig("0.cfg");
  test.ConvertRelativeToAbsolute();
  test.ConvertAbsoluteToRelative();

  // test.UpdateNeighbors(3.5,4.5);
  test.WriteConfig();
  test.WritePOSCAR();


  return 0;
}
