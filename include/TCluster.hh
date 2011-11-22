/*
 * TCluster.hh
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#ifndef TCLUSTER_HH_
#define TCLUSTER_HH_

#include <deque>

#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "TSystem.h"
#include "TObject.h"
using namespace std;
class TCluster :public TObject{
public:
	enum calculationMode_t {maxValue=1, chargeWeighted=2};
	TCluster(){numberOfSeeds=0;numberOfHits=0;seedSigma=0;seedSigma=10;hitSigma=7;isSaturated=false;isGoldenGate=false;};
	TCluster(int eventNumber,int seedSigma=10,int hitSigma=7);
	virtual ~TCluster();
	void addChannel(int channel,Float_t signal,Float_t signalInSigma);
	Float_t getPosition();
	void clear();
	bool isLumpy();
	bool isGoldenGateCluster();
	bool hasSaturatedChannels();
	Float_t getCharge();
	void setPositionCalulation(calculationMode_t mode);
	int size();


private:
	deque< pair<int,Float_t> > cluster;
	int numberOfSeeds;
	int numberOfHits;
	int seedSigma;
	int hitSigma;
	bool isSaturated;
	bool isGoldenGate;
	calculationMode_t mode;
	ClassDef(TCluster,2);
};

#endif /* TCLUSTER_HH_ */
