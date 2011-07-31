/*
 * TDiamondTrack.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TDIAMONDTRACK_HH_
#define TDIAMONDTRACK_HH_
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

#include "TDetectorPlane.hh"

typedef unsigned int uint;

using namespace std;

class TDiamondTrack{

   public:
      TDiamondTrack() {};
      TDiamondTrack(int event, TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3);
      TDiamondTrack(int event, TDetectorPlane Det0, TDetectorPlane Det1, TDetectorPlane Det2, TDetectorPlane Det3, TDetectorPlane Dia);
      ~TDiamondTrack();
      TDetectorPlane const GetD0() { return D0;};
      TDetectorPlane const GetD1() { return D1;};
      TDetectorPlane const GetD2() { return D2;};
      TDetectorPlane const GetD3() { return D3;};
      TDetectorPlane const GetD4() { return D4;};
      TDetectorPlane const GetDia() { return D4;};
      TDetectorPlane const GetD(Int_t detector);
      Float_t GetDetectorHitPosition(Int_t det /*det = 0 to 8*/);
      void SetDetectorHitPosition(Int_t det /*det = 0 to 8*/, Float_t pos);
	int GetEventNumber() {return event_number;};
	void SetEventNumber(int event) {event_number = event;};
	bool FakeTrack;

   protected:
	int event_number;
      TDetectorPlane D0;
      TDetectorPlane D1;
      TDetectorPlane D2;
      TDetectorPlane D3;
      TDetectorPlane D4;
};



#endif /* TDIAMONDTRACK_HH_ */
