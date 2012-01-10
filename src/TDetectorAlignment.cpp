/*
 * TDetectorAlignment.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TDetectorAlignment.hh"
using namespace std;

TDetectorAlignment::TDetectorAlignment(){
	verbosity=5;
   nDetectors = 5;

   for(Int_t i=0; i<5; i++) {
      det_x_offset[i] = 0;
      det_y_offset[i] = 0;
      det_z_offset[i] = 0;
      det_phix_offset[i] = 0;
      det_phiy_offset[i] = 0;
   }
}


//Double_t TDetectorAlignment::GetSiResolution() {
//   Int_t deta[] = {1,0,0,0}, detb[] = {0,1,2,3}, detc[] = {3,3,3,2};
//   Double_t zab, zbc, ref_residual_width, numerator=0, denominator=0;
//   for(int det=0; det<8; det++) {
//      if(det%2) ref_residual_width = det_y_resolution[det/2];
//      else ref_residual_width = det_x_resolution[det/2];
//      zab = TMath::Abs(track_storage[0].GetD(deta[det/2]).GetZ() - track_storage[0].GetD(detb[det/2]).GetZ());
//      zbc = TMath::Abs(track_storage[0].GetD(detb[det/2]).GetZ() - track_storage[0].GetD(detc[det/2]).GetZ());
//      numerator += ref_residual_width*ref_residual_width * (1 + (zab*zab + zbc*zbc)/(zab+zbc)/(zab+zbc));
//      denominator += (1 + (zab*zab + zbc*zbc)/(zab+zbc)/(zab+zbc)) * (1 + (zab*zab + zbc*zbc)/(zab+zbc)/(zab+zbc));
//   }
//
//   return TMath::Sqrt(numerator/denominator);
//}


int TDetectorAlignment::getVerbosity() const
{
    return verbosity;
}

void TDetectorAlignment::AddToPhiXOffset(UInt_t plane, Float_t addPhiXOffset)
{
	if(plane<6){
		vecDetPhiXOffset[plane].push_back(addPhiXOffset);
		det_phix_offset[plane]+=addPhiXOffset;
	}
}
void TDetectorAlignment::AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset)
{
	if(plane<6){
		vecDetPhiYOffset[plane].push_back(addPhiYOffset);
		det_phiy_offset[plane]+=addPhiYOffset;
	}
}

void TDetectorAlignment::AddToXOffset(UInt_t plane, Float_t addXOffset)
{
	if(plane<6){
		vecDetXOffset[plane].push_back(addXOffset);
		det_x_offset[plane]+=addXOffset;
	}
}

void TDetectorAlignment::AddToYOffset(UInt_t plane, Float_t addYOffset)
{
	if(plane<6){
		vecDetYOffset[plane].push_back(addYOffset);
		det_y_offset[plane]+=addYOffset;
	}
}

void TDetectorAlignment::setVerbosity(int verbosity)
{
	if(verbosity!=this->verbosity){
		cout<<"TDetectorAlignment::setVerbosity from "<<getVerbosity()<<" to "<<verbosity<<endl;
		this->verbosity = verbosity;
	}
}


