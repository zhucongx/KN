//
// Created by Zhucong Xi on 1/31/20.
//
#include "Atom.h"
#include "Configuration.h"

int main(int argc, char *argv[]) {
  Configuration test;
  test.readConfig("test");
  test.writeConfig("test1");
  test.cnvPrl2Pst();
  test.cnvPst2Prl();


  test.writeConfig("test2");
  return 1;
}