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
using namespace std;
class TCluster {
public:
	enum calculationMode_t {maxValue=1, chargeWeighted=2};
	TCluster(int eventNumber,int seedSigma,int hitSigma);
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

	//whats up with eta?
private:
	deque< pair<int,Float_t> > cluster;
	int numberOfSeeds;
	int numberOfHits;
	int seedSigma;
	int hitSigma;
	bool isSaturated;
	bool isGoldenGate;
};

#endif /* TCLUSTER_HH_ */
