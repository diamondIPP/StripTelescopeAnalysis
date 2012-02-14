//
//  TEvent.hh
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#ifndef TEVENT_HH_
#define TEVENT_HH_

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

#include "TPlaneProperties.hh"
#include "TPlane.hh"

class TEvent:public TObject {
public:
	TEvent(UInt_t nEvent=0);
	TEvent(const TEvent& rhs);
	virtual ~TEvent();
	TPlane getPlane(int plane){return planes[plane];};
	void addPlane(TPlane plane,Int_t pos=-1);
	UInt_t getNXClusters(UInt_t plane);//
	UInt_t getNYClusters(UInt_t plane);//{if(plane<planes.size())return planes.at(plane).getNYClusters();else return 0;}
	UInt_t getNClusters(UInt_t det);
	TCluster getCluster(UInt_t plane,TPlane::enumCoordinate cor, UInt_t cl);
	TCluster getCluster(UInt_t det, UInt_t cl);
	UInt_t getClusterSize(UInt_t det,UInt_t cl);
	UInt_t getClusterSize(UInt_t plane,TPlane::enumCoordinate cor, UInt_t cl);

	UInt_t getNPlanes();
	bool  isValidSiliconEvent();
	bool isMasked();
	UInt_t getEventNumber(){return eventNumber;};
	void setVerbosity(UInt_t verbosity);
	void Print(UInt_t level);
private:
	
	vector<TPlane> planes;
	UInt_t eventNumber;
	UInt_t verbosity;

    ClassDef(TEvent,6);
public:
	
};

#endif //TEVENT_HH_
