//
// Created by Zhucong Xi on 1/31/20.
//
#include "FCCConfig.h"
int main(int argc, char *argv[]) {
  FCCConfig test;
  // test.GenerateFCC(16.184,"Al",{1,1,1});

  // test.GenerateHCP(3.209,5.211,"Mg",{50,50,50});

  // if(!test.ReadPOSCAR("test.pos")) { return 5;}
  // test.Perturb();
  test.ReadConfig("0.cfg");
  test.UpdateNeighbors(3.5,4.5);
  test.WriteConfig();
  test.WritePOSCAR();


  return 0;
}
