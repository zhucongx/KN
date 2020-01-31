//
// Created by Zhucong Xi on 1/31/20.
//

#ifndef _ATOMCLASS_H_
#define _ATOMCLASS_H_

// For FCC the first nearest neighbor is 12
#define firstNearestNbr 12

class AtomClass {
 private:
  // atom id which is an unique int for every atom
  int id;
  string type;
  // real position
  double position[3];
  // relative position in the box
  double relativePosition[3];
  //
  vector<int> nearNbrList;
  // First nearest neighbor list
  array<int, firstNearestNbr> firstNearestNbrList;
};

#endif //_ATOMCLASS_H_
