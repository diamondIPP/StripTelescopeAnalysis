/*
 * TResidual.hh
 *
 *  Created on: Jan 11, 2012
 *      Author: bachmair
 */

#ifndef TRESIDUAL_HH_
#define TRESIDUAL_HH_

#include "TROOT.h"
#include "TMath.h"
#include "TPlane.hh"
#include "TCluster.hh"
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

class TResidual {
public:
	TResidual();
	virtual ~TResidual();
	void Print(UInt_t level=0);
	void addDataPoint(Float_t deltaX,Float_t predX,Float_t deltaY,Float_t predY);
	Float_t getXMean();
	Float_t getYMean();
	Float_t getXSigma();
	Float_t getYSigma();
	UInt_t getUsedTracks(){return nUsedTracks;}

	Float_t getXOffset();
	Float_t getYOffset();

	Float_t getPhiXOffset();
	Float_t getPhiYOffset();
	void SetTestResidual(bool value=true){bTestResidual=value;}
private:
	Float_t resXMean, resXSigma,resYMean,resYSigma;
	Float_t sumRx;
	Float_t sumRy;
	Float_t sumVx;
	Float_t sumVy;
	Float_t sumV2x;
	Float_t sumV2y;
	Float_t sumVRx;
	Float_t sumVRy;
	UInt_t nUsedTracks;

	bool bTestResidual;
};

#endif /* TRESIDUAL_HH_ */
