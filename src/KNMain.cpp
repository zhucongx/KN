//
// Created by Zhucong Xi on 1/31/20.
//
#include "Atom.h"
#include "Configuration.h"

int main(int argc, char *argv[]) {
  Configuration test;
  test.readConfig("test.conf");
  // test.writeConfig("test1");
  // test.cnvPrl2Pst();
  // test.cnvPst2Prl();

  // test.readPOSCAR("test.pos");
  test.writePOSCAR("test2",false);
  return 1;
}

/***
vacO X
T    T    T
T    F

***/