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
      Double_t GetPhiXOffset(Int_t plane) {return det_phix_offset[plane];};
      void AddToPhiXOffset(UInt_t plane, Float_t addPhiXOffset);
      Double_t GetPhiYOffset(Int_t plane) {return det_phiy_offset[plane];};
      void AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset);//{if(plane<6)det_phiy_offset[plane]+=addPhiYOffset;}
      std::vector<Double_t> GetXOffsetHistory(Int_t plane) {return det_x_offset_history[plane];};
      std::vector<Double_t> GetYOffsetHistory(Int_t plane) {return det_y_offset_history[plane];};
      std::vector<Double_t> GetZOffsetHistory(Int_t plane) {return det_z_offset_history[plane];};
      std::vector<Double_t> GetPhiOffsetHistory(Int_t plane) {return det_phi_offset_history[plane];};
      int getVerbosity() const;
      void setVerbosity(int verbosity);


   private:

      //store global offsets here
      Double_t det_x_offset[6];
      std::vector<Double_t> vecDetXOffset[6];
      Double_t det_y_offset[6];
      std::vector<Double_t> vecDetYOffset[6];
      Double_t det_z_offset[6];
      std::vector<Double_t> vecDetZOffset[6];
      Double_t det_phix_offset[6];
      std::vector<Double_t> vecDetPhiXOffset[6];
      Double_t det_phiy_offset[6];
      std::vector<Double_t> vecDetPhiYOffset[6];

      //store reconstructed offsets here
      std:: vector<Double_t> det_x_offset_history[6];
      std::vector<Double_t> det_y_offset_history[6];
      std::vector<Double_t> det_z_offset_history[6];
      std::vector<Double_t> det_phi_offset_history[6];

      UInt_t nDetectors;
   private:
      int verbosity;


};
#endif /* TDETECTORALIGNMENT_HH_ */
