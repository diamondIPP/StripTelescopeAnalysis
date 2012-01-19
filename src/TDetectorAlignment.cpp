/*
 * TDetectorAlignment.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TDetectorAlignment.hh"
using namespace std;
ClassImp(TDetectorAlignment);
/*
 *   XY		 YX					   XY		YX
 *   ||      ||				       ||       ||
 *   01   	 23					   45	    67
 *   0		  1						2 		3
 *   D0Z	D1Z
 *   0
 *
 *   |||||
 *   -----
 *   PLANE1: X OFF,YOFF,PHIX,PHIY
 */

TDetectorAlignment::TDetectorAlignment(){
	this->SetName("currentAlignment");
	this->SetTitle("currentAlignment");
	verbosity=0;
   nDetectors = 5;

   for(Int_t i=0; i<5; i++) {
      det_x_offset[i] = 0;
      det_y_offset[i] = 0;
      det_z_offset[i] = 0;
      det_phix_offset[i] = 0;
      det_phiy_offset[i] = 0;
      xResolution[i]=0.3;
      yResolution[i]=0.3;
   }
}

//TODO: where to we put get SiResolution...
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
	if(verbosity)cout<<"TDetectorAlignment::addPhiXOffset of Plane "<<plane<<": "<<addPhiXOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetPhiXOffset[plane].push_back(addPhiXOffset);
		det_phix_offset[plane]+=addPhiXOffset;
	}
	if(verbosity)cout<<det_phix_offset[plane]<<endl;
}
void TDetectorAlignment::AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::addPhiYOffset of Plane "<<plane<<": "<<addPhiYOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetPhiYOffset[plane].push_back(addPhiYOffset);
		det_phiy_offset[plane]+=addPhiYOffset;
	}
	if(verbosity)cout<<det_phiy_offset[plane]<<endl;
}

void TDetectorAlignment::AddToXOffset(UInt_t plane, Float_t addXOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::AddToXOffset of Plane "<<plane<<": "<<addXOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetXOffset[plane].push_back(addXOffset);
		det_x_offset[plane]+=addXOffset;
	}
	if(verbosity)cout<<det_x_offset[plane]<<endl;
}

void TDetectorAlignment::AddToYOffset(UInt_t plane, Float_t addYOffset)
{
	if(verbosity)cout<<"TDetectorAlignment::AddToYOffset of Plane "<<plane<<": "<<addYOffset<<", new Value: "<<flush;
	if(plane<6){
		vecDetYOffset[plane].push_back(addYOffset);
		det_y_offset[plane]+=addYOffset;
	}
	if(verbosity)cout<<det_y_offset[plane]<<endl;
}

void TDetectorAlignment::AddToZOffset(UInt_t plane, Float_t addZOffset)
{
	if(plane<6){
		vecDetZOffset[plane].push_back(addZOffset);
		det_z_offset[plane]+=addZOffset;
	}
}

void TDetectorAlignment::PrintXOffset(UInt_t plane,UInt_t level)
{
	cout<<TCluster::Intent(level)<<"X-Offset:" <<det_x_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetXOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		cout<<" "<<vecHis.at(i);
	cout<<endl;
}

void TDetectorAlignment::PrintYOffset(UInt_t plane,UInt_t level)
{
	cout<<TCluster::Intent(level)<<"Y-Offset:" <<det_y_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetYOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		cout<<" "<<vecHis.at(i);
	cout<<endl;
}

void TDetectorAlignment::PrintZOffset(UInt_t plane,UInt_t level)
{
	cout<<TCluster::Intent(level)<<"Z-Offset:" <<det_z_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetZOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		cout<<" "<<vecHis.at(i);
	cout<<endl;
}

void TDetectorAlignment::PrintPhiXOffset(UInt_t plane, UInt_t level)
{
	cout<<TCluster::Intent(level)<<"PhiX-Offset:" <<det_phix_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetPhiXOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		cout<<" "<<vecHis.at(i);
	cout<<endl;
}

void TDetectorAlignment::PrintPhiYOffset(UInt_t plane, UInt_t level)
{
	cout<<TCluster::Intent(level)<<"PhiY-Offset:" <<det_phiy_offset[plane]<<"\t";
	vector<Double_t> vecHis = this->GetPhiYOffsetHistory(plane);
	for(UInt_t i=0;i<vecHis.size();i++)
		cout<<" "<<vecHis.at(i);
	cout<<endl;
}

void TDetectorAlignment::PrintResults(UInt_t level)
{
	cout<<"Alignment Results"<<endl;
	for(UInt_t plane=1;plane<this->nDetectors;plane++){
		cout<<TCluster::Intent(level+1)<<" Plane "<< plane<<endl;
		PrintXOffset(plane,level+2);
		PrintYOffset(plane,level+2);
		PrintZOffset(plane,level+2);
//		cout<<endl;
		PrintPhiXOffset(plane,level+2);
		PrintPhiYOffset(plane,level+2);
		cout<<endl;
	}
}

Double_t TDetectorAlignment::getXResolution(UInt_t plane)
{
	if(plane<6)
    return xResolution[plane];
}

void TDetectorAlignment::setXResolution(Double_t xres,UInt_t plane)
{
	printf("Set X-Resolution of Plane %d to %2.6",plane,xres);
	if(plane<6)
    xResolution[plane] = xres;
}

Double_t TDetectorAlignment::getYResolution(UInt_t plane)
{
	if(plane<6)
    return yResolution[plane];
}

void TDetectorAlignment::setYResolution(Double_t resolution,UInt_t plane)
{
	printf("Set Y-Resolution of Plane %d to %2.6",plane,resolution);
	if(plane<6)
		yResolution[plane] = resolution;
}

void TDetectorAlignment::setVerbosity(int verbosity)
{
	if(verbosity!=this->verbosity){
		cout<<"TDetectorAlignment::setVerbosity from "<<getVerbosity()<<" to "<<verbosity<<endl;
		this->verbosity = verbosity;
	}
}


