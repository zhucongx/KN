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

  test.ReadConfig("test.conf");
  test.ConvertAbsoluteToRelative();
  test.ConvertRelativeToAbsolute();
  test.WriteConfig();
  test.WritePOSCAR();


  return 0;
}
