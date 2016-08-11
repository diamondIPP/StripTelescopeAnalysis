//
//  TPlane.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef TPlane_hh
#define TPlane_hh

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstring>
#include <deque>


//ROOT Class Headers
#include "TROOT.h" // for adding your own classes to ROOT's library
#include "TObject.h"
#include "TH1F.h"

//own classes
#include "TCluster.hh"
#include "TPlaneProperties.hh"


//class TCluster;

class TPlane:public TObject {
public:

	TPlane(){type = TPlaneProperties::kUndefined;xClusters.clear();yClusters.clear();planeNo=0;verbosity=0;};

	TPlane(UInt_t planeNo,vector<TCluster> xClusters, vector<TCluster> yClusters,TPlaneProperties::enumDetectorType type=TPlaneProperties::kSilicon);
	TPlane(UInt_t planeNo,vector<TCluster> xCluster,TPlaneProperties::enumDetectorType type=TPlaneProperties::kDiamond);
//	TPlane(UInt_t planeNo,vector<TCluster> xCluster,TPlaneProperties::enumDetectorType type=TPlaneProperties::kDiamond, Float_t ecmNoise);
//	TPlane(UInt_t planeNo,vector<TCluster> xCluster, TPlaneProperties::enumDetectorType type, Float_t *eadc, Float_t *eped, Float_t *epedCMN, Float_t *epedSigma, Float_t *epedSigmaCMN, Float_t *erawSignal, Float_t *erawSignalCMN, Float_t ecmNoise);
	TPlane(const TPlane& rhs);//COPY Constructor
	TPlane &operator=(const TPlane &src); //class assignment function
	virtual ~TPlane();
	TCluster getCluster(TPlaneProperties::enumCoordinate cor, UInt_t cl);
	TCluster getXCluster(UInt_t cl);
	TCluster getYCluster(UInt_t cl);
	Float_t getXPosition(UInt_t cl,bool cmnCorrected,TCluster::calculationMode_t mode=TCluster::highest2Centroid,TH1F* histo=0);
	Float_t getYPosition(UInt_t cl,bool cmnCorrected,TCluster::calculationMode_t mode=TCluster::highest2Centroid,TH1F* histo=0);
	Float_t getPosition(TPlaneProperties::enumCoordinate cor, UInt_t cl,bool cmnCorrected,TCluster::calculationMode_t mode=TCluster::highest2Centroid,TH1F* histo=0);
	UInt_t getNXClusters();
	UInt_t getNYClusters();
	bool isValidPlane();
    enum TPlaneProperties::enumDetectorType getDetectorType() const;
    void setDetectorType(TPlaneProperties::enumDetectorType type);
    void Print(UInt_t level=0);
    void SetClusters(vector<TCluster> xClusters, vector<TCluster> yClusters);
    void SetXClusters(vector<TCluster> xClusters);
    void SetYClusters(vector<TCluster> yClusters);
    bool hasInvalidReadout();
	void SetSignalValues(UShort_t *eadc, Float_t *eped, Float_t *epedCMN, Float_t *epedSigma, Float_t *epedSigmaCMN, Float_t *erawSignal, Float_t *erawSignalCMN, Float_t ecmNoise);
private:
    TPlaneProperties::enumDetectorType type;
	UInt_t planeNo;
	UInt_t verbosity;
	vector<TCluster> xClusters, yClusters;
	vector<UShort_t> adc;
	vector<Float_t> ped, pedCMN, pedSigma, pedSigmaCMN, rawSignal, rawSignalCMN;
	Float_t cmNoise;
//	vector<Int_t> ch;
    ClassDef(TPlane,9);
};

#endif // TPlane_hh
