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
#include "TROOT.h"
//#define TCLUSTER_REVISION 12;
using namespace std;
class TCluster :public TObject{
public:
	static UInt_t TCLUSTER_REVISION() {return 20;};
    typedef vector<vector<TCluster> > vecvecTCluster;
    enum calculationMode_t{ maxValue = 1, chargeWeighted = 2, highest2Centroid =3};
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
    TCluster(int eventNumber,UChar_t det,  int seedSigma = 10, int hitSigma = 7,UInt_t nChannels=256);
    TCluster(const TCluster& a);
    virtual ~TCluster();
    void addChannel(UInt_t channel, Float_t signal, Float_t signalInSigma, UShort_t adcValue, bool bSaturated,bool isScreened);
    Float_t getPosition(calculationMode_t mode=highest2Centroid);
    void clear();
    bool isLumpyCluster();
    bool isGoldenGateCluster();
    bool hasSaturatedChannels();
    Float_t getCharge(bool useSmallSignals=false);
    Float_t getCharge(UInt_t clusters,bool useSmallSignals=false);
    void setPositionCalulation(calculationMode_t mode);
    UInt_t size();
    UInt_t getHighestSignalChannel();
    Float_t getChargeWeightedMean(bool useNonHits=false);
    void checkCluster();
    bool isSeed(UInt_t cl);
    bool isHit(UInt_t cl);
    Float_t getSignalOfChannel(UInt_t channel);
    UInt_t getSmallestChannelNumber();
    UInt_t getHighestChannelNumber();
    Float_t getHighestSignal();
    Float_t getSignal(UInt_t clusterPos);
    Float_t getSNR(UInt_t clusterPos);
    UShort_t getAdcValue(UInt_t clusterPos);
    UInt_t getHighestHitClusterPosition();
    UInt_t getClusterPosition(UInt_t channelNo);
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
    Float_t getHighest2Centroid();
    void Print(UInt_t level=0);
    static string Intent(UInt_t level);
	Float_t getEta();

private:
    void checkForGoldenGate();
    void checkForLumpyCluster();
    UInt_t checkClusterForSize() const;
    //deque<pair<int,Float_t> > cluster; //ch,signal
    deque <UInt_t> clusterChannel;
    deque <Float_t> clusterSignal;
    //        deque<pair<UShort_t,Float_t> > cluster2; //adc,SNR
    deque<UShort_t> clusterADC;
    deque<Float_t> clusterSignalInSigma;
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
    UChar_t det;
    UInt_t eventNumber;
    ClassDef(TCluster,TCLUSTER_REVISION());
};

#endif /* TCLUSTER_HH_ */
