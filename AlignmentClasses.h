//Classes used for silicon and diamond alignment
//based on Alignment_Monte_Carlo_j15a.cpp
//crapload of changes to class structures
//fixed segafaulting in phi and z alignment by moving plot generation to its own member fcn
//2010-09-18 Plucked from dev-2010-04-30-eta_correction for adapatation in dev-2010-09-18-alignment
//2010-09
//           TODO: Several improvements would be tracking and aligning x and y planes separately and of course tearing the damn telescope apart to pin down the geometry
//New comment
#ifndef ALIGNMENTCLASSES_H
#define ALIGNMENTCLASSES_H
//C++ Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
//#include <deque>
#include <ctime> // seed the random number generator
#include <cstdlib> // random number generator
#include <sstream>

#include "TDiamondTrack.hh"
#include "TDetectorPlane.hh"
#include "TDetectorAlignment.hh"
#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

typedef unsigned int uint;

using namespace std;
//using namespace TMath;


#endif /*ALIGNMENTCLASSES_H*/
