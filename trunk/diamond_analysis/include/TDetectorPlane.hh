/*
 * TDetectorPlane.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TDETEKTORPLANE_HH_
#define TDETEKTORPLANE_HH_
//C++ Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <ctime> // seed the random number generator
#include <cstdlib> // random number generator
#include <sstream>

//Root Libraries
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
class TDetectorPlane {

   public:
      TDetectorPlane() {};
      ~TDetectorPlane() {};
      void SetX(Float_t x) {X_position = x;};
      void SetY(Float_t y) {Y_position = y;};
      void SetZ(Float_t z) {Z_position = z;};
      Float_t GetX() const {return X_position;}
      Float_t GetY() const {return Y_position;}
      Float_t GetZ() const {return Z_position;}

   private:
      Float_t X_position;
      Float_t Y_position;
      Float_t Z_position;

};

#endif /* TDETEKTORPLANE_HH_ */
