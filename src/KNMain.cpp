//
// Created by Zhucong Xi on 1/31/20.
//
#include "Config.h"

int main(int argc, char *argv[]) {
  Config test;
  // test.generateFCC(16.184,"Al",std::vector<int>(3,10));
  // test.readPOSCAR("test.pos");
  // test.cnvPrl2Pst();
  // test.cnvPst2Prl();

  // if(!test.ReadPOSCAR("test.pos")) { return 5;}
  // test.ConvertAbsoluteToRelative();
  // test.ConvertRelativeToAbsolute();
  test.GenerateHCP(3.209,5.211,"Mg",{2,2,2});
  test.Perturb();
  test.WriteConfig();
  test.WritePOSCAR();


  return 0;
}
