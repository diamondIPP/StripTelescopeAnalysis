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

//own classes
#include "TCluster.hh"


#define N_INVALID -9999
//class TCluster;
class TPlane:public TObject {
public:

	TPlane(){type = kUndefined;};
	enum enumDetectorType{kUndefined = 0, kSilicon = 1, kDiamond =2};
	enum enumCoordinate{ X_COR =1, Y_COR=2, Z_COR =3, XY_COR=4,};
	TPlane(vector<TCluster> xClusters, vector<TCluster> yClusters,enumDetectorType type=kSilicon);
	TPlane(vector<TCluster> xCluster,enumDetectorType type=kDiamond);
	virtual ~TPlane();
	Float_t getXPosition(UInt_t cl,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getYPosition(UInt_t cl,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getPosition(enumCoordinate cor, UInt_t cl,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	UInt_t getNXClusters();
	UInt_t getNYClusters();
	bool isValidPlane();
    enum enumDetectorType getDetectorType() const;
    void setDetectorType(enumDetectorType type);
    static string getCoordinateString(enumCoordinate cor);
    static string getDetectortypeString(enumDetectorType type);
    void Print(UInt_t level=0);
private:
	enumDetectorType type;
	
	vector<TCluster> xClusters, yClusters;
    ClassDef(TPlane,1);
};

#endif // TPlane_hh
