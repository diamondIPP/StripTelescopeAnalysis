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
#include "TPlane.hh"
//
//#include "TRandom3.h"
//#include "TMath.h"
//#include "TF1.h"
//#include "TGraph.h"
//#include "TCanvas.h"
//#include "TH1F.h"
#include "TROOT.h"
#include "TDatime.h"
#include "TNamed.h"
//#include "TApplication.h"
//#include "TSystem.h"
//#include "TH2F.h"
//#include "TFile.h"
//#include "TStyle.h"
//
//#include "TDiamondTrack.hh"
//#include "TDetectorPlane.hh"
#include "TCluster.hh"


class TDetectorAlignment:public TNamed{

public:
  TDetectorAlignment();
  ~TDetectorAlignment() {};
  TDetectorAlignment &operator=(const TDetectorAlignment &src);
  void Print(Option_t *opt = "");
  Double_t GetXOffset(UInt_t plane) {return det_x_offset[plane];};
  Double_t GetYOffset(UInt_t plane) {return det_y_offset[plane];};
  Double_t GetZOffset(UInt_t plane) {return det_z_offset[plane];};

  void SetXOffset(UInt_t plane,Float_t xOffset) {if(plane<nDetectors)det_x_offset[plane]=xOffset;};
  void SetYOffset(UInt_t plane,Float_t yOffset) {if(plane<nDetectors)det_y_offset[plane]=yOffset;};
  void SetZOffset(UInt_t plane,Float_t zOffset) {if(plane<nDetectors)det_z_offset[plane]=zOffset;};

  void AddToXOffset(UInt_t plane, Float_t addXOffset);//{if(plane<6)det_x_offset[plane]+=addXOffset;}
  void AddToYOffset(UInt_t plane, Float_t addYOffset);//{if(plane<6)det_y_offset[plane]+=addYOffset;}
  void AddToZOffset(UInt_t plane, Float_t addZOffset);

  Double_t GetPhiXOffset(UInt_t plane) {return det_phix_offset[plane];};
  Double_t GetPhiYOffset(UInt_t plane) {return det_phiy_offset[plane];};

  void AddToPhiXOffset(UInt_t plane, Float_t addPhiXOffset);
  void AddToPhiYOffset(UInt_t plane, Float_t addPhiYOffset);

  Double_t GetLastXOffset(UInt_t plane)  {if (plane<nDetectors)return vecDetXOffset[plane].size()==0?0:vecDetXOffset[plane].back();return -9999;};
  Double_t GetLastYOffset(UInt_t plane)  {if (plane<nDetectors)return vecDetYOffset[plane].size()==0?0:vecDetYOffset[plane].back();return -9999;};
  Double_t GetLastPhiXOffset(UInt_t plane) {if (plane<nDetectors)return vecDetPhiXOffset[plane].size()==0?0:vecDetPhiXOffset[plane].back();return -9999;};
  Double_t GetLastPhiYOffset(UInt_t plane) {if (plane<nDetectors)return vecDetPhiYOffset[plane].size()==0?0:vecDetPhiXOffset[plane].back();return -9999;};

  std::vector<Double_t> GetXOffsetHistory(UInt_t plane) {if (plane<nDetectors)return vecDetXOffset[plane];std::vector<Double_t> a;return a;};
  std::vector<Double_t> GetYOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetYOffset[plane];std::vector<Double_t> a;return a;};
  std::vector<Double_t> GetZOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetZOffset[plane];std::vector<Double_t> a;return a;};
  std::vector<Double_t> GetPhiXOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetPhiXOffset[plane];std::vector<Double_t> a;return a;};
  std::vector<Double_t> GetPhiYOffsetHistory(UInt_t plane) {if (plane<nDetectors) return vecDetPhiYOffset[plane];std::vector<Double_t> a;return a;};

  std::string PrintXOffset(UInt_t plane,UInt_t level);
  std::string PrintYOffset(UInt_t plane,UInt_t level);
  std::string PrintZOffset(UInt_t plane,UInt_t level);

  std::string PrintPhiXOffset(UInt_t plane,UInt_t level);
  std::string PrintPhiYOffset(UInt_t plane,UInt_t level);

  std::string PrintResolution(UInt_t plane,UInt_t level);


  std::string PrintResults(UInt_t level=0);
  void ResetAlignment(Int_t plane =-1);
  int getVerbosity() const;
  void setVerbosity(int verbosity);
  Double_t getXResolution(UInt_t plane);
  Double_t getZResolution(UInt_t plane);
  void setXResolution(Double_t xRes,UInt_t plane);
  void setResolution(Double_t xRes, UInt_t plane,TPlaneProperties::enumCoordinate cor);
  Double_t getXMean(UInt_t plane);
  void setXMean(Double_t xMean,UInt_t plane);

  Double_t getYResolution(UInt_t plane);
  void setYResolution(Double_t yRes,UInt_t plane);
  Double_t getResolution(UInt_t plane, TPlaneProperties::enumCoordinate cor){return cor==TPlaneProperties::X_COR?getXResolution(plane):getXResolution(plane);}
  Double_t getYMean(UInt_t plane);
  void setYMean(Double_t yMean,UInt_t plane);
  UInt_t getNUsedEvents() const {return (this->nUsedEvents);};
  void setNUsedEvents(UInt_t usedEvents){this->nUsedEvents=usedEvents;};
  bool isPreAligned(Float_t maxOffset=0.3,Int_t nAlignedDetectors=2);
  void addEventIntervall(UInt_t first,UInt_t last);
  void setAlignmentTrainingTrackFraction(Float_t fraction){alignmentTrainingTrackFraction=fraction;};
  void setDiamondDate(){diaTime=TDatime();};
  void setSiliconDate(){silTime=TDatime();};
  void setRunNumber(UInt_t rn){this->runNumber=rn;};
  void setDiaChi2(Float_t chi2){this->diaChi2=chi2;};
  Float_t getDiaChi2(){return diaChi2;};
  void setNDiamondAlignmentEvents(UInt_t nEvents){this->nDiamondAlignmentEvents=nEvents;};
  UInt_t getDiamondAlignmentEvents(){return nDiamondAlignmentEvents;};
  string getLastUpdateTimeAsString(){
    if(diaTime>silTime)return diaTime.AsString();
    else return silTime.AsString();

  }
private:
  void UpdateTime(UInt_t plane){if(plane<4)setSiliconDate();else setDiamondDate();};
  Double_t xResolution[6];//Planes
  Double_t yResolution[6];
  Double_t xMean[6];//Planes
    Double_t yMean[6];
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
  UInt_t runNumber;
  TDatime diaTime,silTime;
  UInt_t nUsedEvents;
  Float_t alignmentTrainingTrackFraction;
  Float_t diaChi2;
  std::vector< UInt_t> intervallBeginEventNo;
  std::vector<UInt_t > intervallEndEventNo;
  UInt_t nDiamondAlignmentEvents;
  UInt_t nEvents;
private:
  int verbosity;
  ClassDef(TDetectorAlignment,2);
};
#endif /* TDETECTORALIGNMENT_HH_ */
