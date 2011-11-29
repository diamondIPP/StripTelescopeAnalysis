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
//#define TCLUSTER_REVISION 7;
using namespace std;
class TCluster :public TObject{
public:
	static UInt_t TCLUSTER_REVISION() {return 9;};
    typedef vector<vector<TCluster> > vecvecTCluster;
    enum calculationMode_t{ maxValue = 1, chargeWeighted = 2};
    TCluster()
    {
        numberOfSeeds = 0;
        numberOfHits = 0;
        seedSigma = 0;
        seedSigma = 10;
        hitSigma = 7;
        isSaturated = false;
        isGoldenGate = false;
        isLumpy = false;
        verbosity = 0;
        maximumSignal = 0;
        charge = 0;
        revisionNumber=TCLUSTER_REVISION();
        isChecked = false;
        hasBadChannel=false;
        numberOfNoHits=0;
        nChannels=256;
    };
    TCluster(int eventNumber, int seedSigma = 10, int hitSigma = 7,UInt_t nChannels=256);
    virtual ~TCluster();
    void addChannel(int channel, Float_t signal, Float_t signalInSigma, UShort_t adcValue, bool bSaturated,bool isScreened);
    Float_t getPosition();
    void clear();
    bool isLumpyCluster();
    bool isGoldenGateCluster();
    bool hasSaturatedChannels();
    Float_t getCharge();
    void setPositionCalulation(calculationMode_t mode);
    UInt_t size();
    int getMaximumChannel();
    Float_t getChargeWeightedMean();
    void checkCluster();
    bool isSeed(UInt_t cl);
    UInt_t getMinChannelNumber();
    UInt_t getMaxChannelNumber();
    Float_t getSignal(UInt_t clusterPos);
    Float_t getSNR(UInt_t clusterPos);
    UShort_t getAdcValue(UInt_t clusterPos);
    UInt_t getChannel(UInt_t clusterPos);
    Float_t getPedestalSigma(UInt_t clusterPos);
    Float_t getPedestalMean(UInt_t clusterPos);
    int getHitSigma() const;
    int getSeedSigma() const;
    void setHitSigma(int hitSigma);
    void setSeedSigma(int seedSigma);
    bool isSaturatedCluster(){return isSaturated;};
    bool isBadChannelCluster(){return hasBadChannel;}
    bool isScreened();
    bool isScreened(UInt_t cl);
    Float_t highest2_centroid();

private:
    void checkForGoldenGate();
    void checkForLumpyCluster();
    deque<pair<int,Float_t> > cluster; //ch,signal
    deque<pair<UShort_t,Float_t> > cluster2; //adc,SNR
    deque<bool> clusterChannelScreened;
    UInt_t numberOfSeeds;
    UInt_t numberOfHits;
    UInt_t numberOfNoHits;
    int seedSigma;
    int hitSigma;
    bool isSaturated;
    bool isGoldenGate;
    bool isLumpy;
    bool isChecked;
    bool hasBadChannel;
    calculationMode_t mode;
    int verbosity;
    Float_t charge;
    Float_t maximumSignal;
    int maxChannel;
    int revisionNumber;
    UInt_t nChannels;
    ClassDef(TCluster,TCLUSTER_REVISION());
};

#endif /* TCLUSTER_HH_ */
