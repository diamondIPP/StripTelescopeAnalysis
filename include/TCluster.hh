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
#define TCLUSTER_REVISION 6 //UPDATE REVISION IF YOU MAKE CHANGES IN TCLUSTER!!!!!
public:
	typedef  vector< vector<TCluster> > vecvecTCluster;
	enum calculationMode_t {maxValue=1, chargeWeighted=2};
	TCluster(){numberOfSeeds=0;numberOfHits=0;seedSigma=0;seedSigma=10;hitSigma=7;isSaturated=false;isGoldenGate=false;isLumpy=false;verbosity=0;maximumSignal=0;charge=0;revisionNumber=TCLUSTER_REVISION;isChecked=false;};
	TCluster(int eventNumber,int seedSigma=10,int hitSigma=7);
	virtual ~TCluster();
	void addChannel(int channel,Float_t signal,Float_t signalInSigma,UShort_t adcValue,bool bSaturated);
	Float_t getPosition();
	void clear();
	bool isLumpyCluster();
	bool isGoldenGateCluster();
	bool hasSaturatedChannels();
	Float_t getCharge();
	void setPositionCalulation(calculationMode_t mode);
	int size();
	int getMaximumChannel();
	Float_t getChargeWeightedMean();
	void checkCluster();
	bool isSeed(UInt_t cl);
	UInt_t getMinChannelNumber();
	UInt_t getMaxChannelNumber();


private:
	void checkForGoldenGate();
	void checkForLumpyCluster();
	deque< pair< int, Float_t> > cluster;
	deque< pair< UShort_t, Float_t> > cluster2;
	int numberOfSeeds;
	int numberOfHits;
	int seedSigma;
	int hitSigma;
	bool isSaturated;
	bool isGoldenGate;
	bool isLumpy;
	bool isChecked;
	calculationMode_t mode;
	int verbosity;
	Float_t charge;
	Float_t  maximumSignal;
	int maxChannel;
	int revisionNumber;
	ClassDef(TCluster,TCLUSTER_REVISION);
};

#endif /* TCLUSTER_HH_ */
