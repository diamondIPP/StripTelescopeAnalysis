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
/** Residual Calculation Class
 *
 * This class is used to calculated the values to correct the current Alignment
 * With adding datapoint the residual will be updated
 * The current correction results can be read out all the time
 * @author Felix Bachmair
 * @date	12.1.2012	last update
 * @date	10.1.2012	creation
 *
 */
class TResidual {
public:
	TResidual(bool bTest=false);
	virtual ~TResidual();
	void Print(UInt_t level=0);
	void addDataPoint(Float_t deltaX,Float_t predX,Float_t deltaY,Float_t predY);
	Float_t getXMean();
	Float_t getYMean();
	Float_t getXSigma();
	Float_t getYSigma();

	Float_t getDeltaXMean();
	Float_t getDeltaYMean();

	Float_t getPredXMean();
	Float_t getPredYMean();

	UInt_t getUsedTracks(){return nUsedTracks;}

	Float_t getXOffset();
	Float_t getYOffset();

//
//	Float_t getLinRegOffsetX();
//	Float_t getLinRegOffsetY();
//	Float_t getLinRegSlopeX();
//	Float_t getLinRegSlopeY();


	Float_t getPhiXOffset();
	Float_t getPhiYOffset();
	void SetTestResidual(bool value=true){bTestResidual=value;}

private:
	std::vector<Float_t> vecDeltaX;
	std::vector<Float_t> vecPredX;
	std::vector<Float_t> vecDeltaY;
	std::vector<Float_t> vecPredY;
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
