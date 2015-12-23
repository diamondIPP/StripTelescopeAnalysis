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
#include "TPlaneProperties.hh"
#include "TH1F.h"
#include <math.h>  
#define TCLUSTER_REV  40
#define INVALID_POSITION -9999
using namespace std;
class TCluster :public TObject{
public:
    static UInt_t TCLUSTER_REVISION() {return TCLUSTER_REV;};
    enum calculationMode_t{ maxValue = 1, chargeWeighted = 2, highest2Centroid = 3,eta=4,corEta=5, highest2CentroidNoSmallHits=6};
    enum direction_t {left = -1,right = +1};
    TCluster()
    {	initialiseNewCluster();
    };
    TCluster(int eventNumber,UChar_t det,  int seedSigma = 10, int hitSigma = 7,UInt_t nChannels=256, float cmNoise=0);
    TCluster(const TCluster& a);//COPY Constructor
    virtual ~TCluster();
    TCluster &operator=(const TCluster &src); //class assignment function

    void addChannel(UInt_t channel, Float_t pedMean, Float_t pedSigma, Float_t pedMeanCMN, Float_t pedSigmaCMN, Int_t adcValue, bool bSaturated,bool isScreened);
    //    void addChannel(UInt_t channel, Float_t signal, Float_t signalInSigma, UShort_t adcValue, bool bSaturated,bool isScreened);
    Float_t getPosition(bool cmnCorrected, calculationMode_t mode=highest2Centroid,TH1F *histo=0);
    void clear();
    bool isLumpyCluster();
    bool isGoldenGateCluster();
    bool hasSaturatedChannels();
    UInt_t getDetector(){return this->det;}
    TCluster getCrossTalkCorrectedCluster(Float_t alpha);

    Float_t getCharge(bool cmnCorrected = false,bool useSmallSignals=false);
    Float_t getCharge(UInt_t clusters,bool cmnCorrected = false,bool useSmallSignals=false);
    Float_t getCharge(int clusters,bool cmnCorrected = false,bool useSmallSignals=false){return getCharge((UInt_t)clusters,cmnCorrected,useSmallSignals);};


    Float_t getPositiveCharge(bool cmnCorrected = false,bool useSmallSignals=false);
    Float_t getPositiveCharge(int clusters,bool cmnCorrected = false,bool useSmallSignals=false){return getPositiveCharge((UInt_t)clusters,cmnCorrected,useSmallSignals);};
    Float_t getPositiveCharge(UInt_t clusters,bool cmnCorrected = false,bool useSmallSignals=false);

    Float_t getTransparentCharge(UInt_t clusters,bool cmnCorrected, bool useSmallSignals);

    bool hasNegativeCharge(Float_t &charge, Int_t &clPos,bool cmnCorrected = false);
    void setPositionCalulation(calculationMode_t mode);
    UInt_t size();
    UInt_t seedSize();
    void UpdateHighestSignalChannel();
    UInt_t getHighestSignalChannel();
    UInt_t GetHighestSignalChannelTransparentCluster();
    UInt_t getHighestSignalNeighbourChannel(UInt_t channelNo,bool cmnCorrected=false);
    UInt_t getHighestSignalNeighbourClusterPosition(UInt_t clPos,bool cmnCorrected=false,bool bNegativeSignals = false);
    Float_t getChargeWeightedMean(bool cmnCorrected, bool useNonHits=false);
    Float_t getEtaPostion(bool cmnCorrected=false);
    Float_t getPositionCorEta(bool cmnCorrected,TH1F* histo=0);
    void checkCluster();
    bool isSeed(UInt_t cl);
    bool isHit(UInt_t cl);
    Float_t getSignalOfChannel(UInt_t channe, bool cmnCorrected=false);
    UInt_t getSmallestChannelNumber();
    UInt_t getHighestChannelNumber();
    Float_t getHighestSignal(bool cmnCorrected=false);
    Float_t getSignal(UInt_t clusterPos, bool cmnCorrected=false);
    Float_t getSNR(UInt_t clusterPos, bool cmnCorrected=false);
    Float_t getCMN(){return cmNoise;}
    Float_t getPedestalMean(UInt_t clusterPos, bool cmnCorrected=false);
    Float_t getPedestalSigma(UInt_t clusterPos, bool cmnCorrected=false);
    Int_t getAdcValue(UInt_t clusterPos);
    UInt_t getHighestHitClusterPosition();
    UInt_t getHighestHitClusterPositionTransparentCluster();
    UInt_t getClusterPosition(UInt_t channelNo);
    UInt_t getChannel(UInt_t clusterPos);
    UInt_t getFirstHitChannel();
    UInt_t getLastHitChannel();
    //    Float_t getPedestalSigma(UInt_t clusterPos);
    //    Float_t getPedestalMean(UInt_t clusterPos);
    int getHitSigma() const;
    int getSeedSigma() const;
    void setHitSigma(int hitSigma);
    void setSeedSigma(int seedSigma);
    bool isSaturatedCluster(){return isSaturated;};
    bool isBadChannelCluster(){return hasBadChannel;}
    bool isScreened();
    bool isScreened(UInt_t cl);
    Float_t getHighest2Centroid(bool cmnCorrected = false, bool useSmallSignals=true);
    void Print(UInt_t level=0);
    static string Intent(UInt_t level);
    Float_t getReversedEta(bool cmnCorrected=false);
    Float_t getReversedEta(Int_t &rightChannel,bool cmnCorrected=false);
    Float_t getEta(bool cmnCorrected=false);
    Float_t getEta(Int_t &leftChannel,bool cmnCorrected=false);
    Float_t getEta(UInt_t clusPos,Int_t &leftChannel,bool cmnCorrected=false);
    UInt_t getClusterSize();
    void setVerbosity(UInt_t verbosity){this->verbosity=verbosity;};
    Float_t getLeftEta(bool cmnCorrected=false);
    Float_t getRightEta(bool cmnCorrected=false);
    static Float_t getValueOfHisto(Float_t x, TH1F* histo);
    UInt_t getEventNumber(){return eventNumber;};
    bool hasInvalidReadout();
    bool IsTransparentCluster(){return !(isTransparentCluster<0);}
    Float_t GetTransparentHitPosition(){return isTransparentCluster;}
    void SetTransparentCluster(Float_t startChannel);
    void SetTransparentClusterSize(UInt_t size);
    UInt_t GetTransparentClusterSize(){return transparentClusterSize;}
    Int_t getTransparentClusterPosition(UInt_t clusterNo=0);
    bool IsValidTransparentClusterPosition(UInt_t clusterPosition);
private:
    direction_t getDirection(UInt_t clusterPosition);
    Float_t getChargeStartingAt(UInt_t nChannels,UInt_t startingClusterPos,direction_t direction, bool noNegativeCharge, bool useCMcorrection, bool useSmallSignals);
    void initialiseNewCluster();
    void checkForGoldenGate();
    void checkForLumpyCluster();
    UInt_t checkClusterForSize() const;
    deque <UInt_t> clusterChannel;
    deque <Float_t> clusterPedMean;
    deque <Float_t> clusterPedMeanCMN;
    deque <Float_t> clusterPedSigma;
    deque <Float_t> clusterPedSigmaCMN;
    deque<Int_t> clusterADC;

    deque <Float_t> clusterSignal;
    deque <Float_t> clusterSignalCMN;
    deque<Float_t> clusterSignalInSigma;
    deque<Float_t> clusterSignalInSigmaCMN;
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
    Float_t cmNoise;
    Float_t isTransparentCluster;
    UInt_t transparentClusterSize;
    ClassDef(TCluster,TCLUSTER_REV);
};
#endif /* TCLUSTER_HH_ */
