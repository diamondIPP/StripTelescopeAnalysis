/*
 * TDetectorAlignment.hh
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#ifndef TDETECTORALIGNMENT_HH_
#define TDETECTORALIGNMENT_HH_
//C++ Libraries
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <iomanip>
//#include <deque>
#include <ctime> // seed the random number generator
#include <cstdlib> // random number generator
#include <sstream>

#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"

#include "TDiamondTrack.hh"
#include "TDetectorPlane.hh"


class TDetectorAlignment{

   public:
	TDetectorAlignment();
      ~TDetectorAlignment() {};

      Double_t GetXOffset(Int_t plane) {return det_x_offset[plane];};
      Double_t GetYOffset(Int_t plane) {return det_y_offset[plane];};
      Double_t GetZOffset(Int_t plane) {return det_z_offset[plane];};

      void SetZOffset(Int_t plane,Float_t zOffset) {if(plane<6)det_z_offset[plane]=zOffset;};

      void AddToXOffset(UInt_t plane, Float_t addXOffset);//{if(plane<6)det_x_offset[plane]+=addXOffset;}
      void AddToYOffset(UInt_t plane, Float_t addYOffset);//{if(plane<6)det_y_offset[plane]+=addYOffset;}
      void AddToZOffset(UInt_t plane, Float_t addZOffset);

      Double_t GetPhiXOffset(Int_t plane) {return det_phix_offset[plane];};
      Double_t GetPhiYOffset(Int_t plane) {return det_phiy_offset[plane];};

      void AddToPhiXOffset(UInt_t plane, Float_t addPhiXOffset);
      void AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset);

      std::vector<Double_t> GetXOffsetHistory(UInt_t plane) {if (plane<nDetectors)return vecDetXOffset[plane];std::vector<Double_t> a;return a;};
      std::vector<Double_t> GetYOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetYOffset[plane];std::vector<Double_t> a;return a;};
      std::vector<Double_t> GetZOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetZOffset[plane];std::vector<Double_t> a;return a;};
      std::vector<Double_t> GetPhiXOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetPhiXOffset[plane];std::vector<Double_t> a;return a;};
      std::vector<Double_t> GetPhiYOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetPhiYOffset[plane];std::vector<Double_t> a;return a;};

      int getVerbosity() const;
      void setVerbosity(int verbosity);


   private:

      //store global offsets here
      Double_t det_x_offset[6];
      Double_t det_y_offset[6];
      Double_t det_z_offset[6];

      std::vector<Double_t> vecDetXOffset[6];
      std::vector<Double_t> vecDetYOffset[6];
      std::vector<Double_t> vecDetZOffset[6];

      Double_t det_phix_offset[6];
      Double_t det_phiy_offset[6];

      std::vector<Double_t> vecDetPhiXOffset[6];
      std::vector<Double_t> vecDetPhiYOffset[6];

      UInt_t nDetectors;
   private:
      int verbosity;


};
#endif /* TDETECTORALIGNMENT_HH_ */
