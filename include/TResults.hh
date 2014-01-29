/*
 * TResults.hh
 *
 *  Created on: May 29, 2012
 *      Author: bachmair
 */

#ifndef TRESULTS_HH_
#define TRESULTS_HH_
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <iostream>
#include <string>
#include <set>
#include <cstring>
//#include "TIter.h"
#include "TKey.h"
#include <TObject.h>
#include <vector>
#include "TPlaneProperties.hh"
#include "TSettings.class.hh"
#include "TDetectorAlignment.hh"
#include <TDatime.h>


class TResults: public TNamed {
public:
    TResults(UInt_t runnumber=0);
    TResults(TSettings *settings);
    TResults(const   TResults& rhs);//copy constructor
    //  TResults &operator=(const   TResults &src); //class assignment function
    virtual ~TResults();
    void createResultFiles(){createOutputTextFile();createOutputResultFile();}
    void saveResults(TString name);
    void openResults(TSettings *Settings);
    void Print();

public:
    void setAlignment(TDetectorAlignment* newAlignment);
    TDatime getLastUpdateDate(){return lastUpdate;};
    void setNoise(UInt_t det, Float_t detNoise);
    void setDiaNoiseCMcorrected(Float_t noise){diaCMCNoise = noise;}
    void setCMN(Float_t cmn);
    UInt_t getAllEvents() const{return nAllEvents;}
    UInt_t getExactlyOneDiamondHit() const { return nExactlyOneDiamondHit;}
    UInt_t getMoreThanOneDiamondHit() const {return nMoreThanOneDiamondHit;}
    UInt_t getNoDiamondHit() const {return nNoDiamondHit;}
    UInt_t getNoSiliconHit() const {return nNoSiliconHit;}
    UInt_t getOneAndOnlyOneSiliconNotFiducialCut() const {return nOneAndOnlyOneSiliconNotFiducialCut;}
    UInt_t getValidSiliconTrack() const {return nValidSiliconTrack;}
    UInt_t getUseForAlignment() const {return nUseForAlignment;}
    UInt_t getUseForAnalysis() const {return nUseForAnalysis;}
    std::pair<Float_t, Float_t> getAvergSiliconCorrection();
    Float_t getAvergDiamondCorrection();
    void setAllEvents(UInt_t allEvents){nAllEvents = allEvents;}
    void setExactlyOneDiamondHit(UInt_t exactlyOneDiamondHit){ nExactlyOneDiamondHit = exactlyOneDiamondHit;}
    void setMoreThanOneDiamondHit(UInt_t moreThanOneDiamondHit){ nMoreThanOneDiamondHit = moreThanOneDiamondHit;}
    void setNoDiamondHit(UInt_t noDiamondHit){ nNoDiamondHit = noDiamondHit;}
    void setNoSiliconHit(UInt_t noSiliconHit){ nNoSiliconHit = noSiliconHit;}
    void setOneAndOnlyOneSiliconNotFiducialCut(UInt_t oneAndOnlyOneSiliconNotFiducialCut){ nOneAndOnlyOneSiliconNotFiducialCut = oneAndOnlyOneSiliconNotFiducialCut;}
    void setValidSiliconTrack(UInt_t validSiliconTrack){ nValidSiliconTrack = validSiliconTrack;};
    void setUseForAlignment(UInt_t useForAlignment){ nUseForAlignment = useForAlignment;}
    void setUseForAnalysis(UInt_t useForAnalysis){nUseForAnalysis = useForAnalysis;};
    void setSeedSigma(UInt_t det,Float_t sigma){if(det<seedSigma.size())seedSigma[det]=sigma;}
    void setHitSigma(UInt_t det,Float_t sigma){if(det<hitSigma.size())hitSigma[det]=sigma;}
    int getRunNumber() const;
    void setRunNumber(UInt_t runNumber);
    void setDoubleGaussianResolution(Float_t gaus1,Float_t gaus2,TSettings::alignmentMode mode);
    void setSingleGaussianResolution(Float_t gaus,TSettings::alignmentMode mode);
    void setSingleGaussianShortResolution(Float_t gaus,TSettings::alignmentMode mode);
    void setSingleGaussianFWTMResolution(Float_t gaus,TSettings::alignmentMode mode);
    void setSingleGaussianFixedResolution(Float_t gaus,TSettings::alignmentMode mode);
    void setSignalFeedOverCorrection(UInt_t det, Float_t correction);
    void setPH_clustered(Float_t mean,Float_t mp,Float_t width,Float_t gSigma, TSettings::alignmentMode mode);
    void updated(){lastUpdate=TDatime();}
public:
    TString getMaskedChannels();
    void initialiseResults();
    void inheritOldResults(const TResults & rhs);
    void setPH_2outOf10(Float_t Mean, Float_t MP, Float_t width, Float_t gSigma,TSettings::alignmentMode mode);
    void setPH_NoutOfN(vector<Float_t> vecMeans, TSettings::alignmentMode mode);
    void createOutputTextFile();
    void createOutputResultFile();
    std::pair<Float_t,Float_t> getAvrgSilNoise();
    string createSection(TString sectionName,map<TString,TString> results);
    void setResultsFromSettings(TSettings * settings);
private:
    string emptyString(UInt_t nChars,char character=' ');
    TString rootFileName;
    TString runDescription;
    TDetectorAlignment alignment;
    //  TSettings Settings;
    TDatime lastUpdate;
    std::string path;
    UInt_t runnumber;
    //  TSettings *settings;
    std::vector<Float_t> seedSigma;
    std::vector<Float_t> hitSigma;
    std::vector<Float_t> noise;
    Float_t diaCMCNoise;
    Float_t CMN;
    std::vector<Float_t> clusterSize;
    std::vector<Float_t> clusterSeedSize;
    std::vector<Float_t> nClusters;
    std::vector<Float_t> clusterPHmean;
    std::vector<Float_t> clusterPHwidth;
    std::vector<Float_t> signalFeedOverCorrection;
    UInt_t nAllEvents;
    UInt_t nNoSiliconHit;
    UInt_t nOneAndOnlyOneSiliconNotFiducialCut;
    UInt_t nValidSiliconTrack;
    UInt_t nNoDiamondHit;
    UInt_t nMoreThanOneDiamondHit;
    UInt_t nExactlyOneDiamondHit;
    UInt_t nUseForAlignment;
    UInt_t nUseForAnalysis;

    Float_t mean_clustered_normal;
    Float_t mp_clustered_normal;
    Float_t width_clustered_normal;
    Float_t gSigma_clustered_normal;

    Float_t mean_clustered_trans;
    Float_t mp_clustered_trans;
    Float_t width_clustered_trans;
    Float_t gSigma_clustered_trans;

    Float_t mean2outOf10_normal;
    Float_t mp2outOf10_normal;
    Float_t width2outOf10_normal;
    Float_t gSigma2outOf10_normal;
    vector<Float_t> meanNoutOfN_normal;

    Float_t mean2outOf10_trans;
    Float_t mp2outOf10_trans;
    Float_t width2outOf10_trans;
    Float_t gSigma2outOf10_trans;
    vector<Float_t> meanNoutOfN_trans;
    TString textFileName;
    TString resultsFileName;

    //Resolutions
    Float_t doubleGaus1_normal;
    Float_t doubleGaus2_normal;
    Float_t doubleGaus1_trans;
    Float_t doubleGaus2_trans;

    Float_t singleGausFixed_normal;
    Float_t singleGausFixed_trans;
    Float_t singleGausShort_normal;
    Float_t singleGausShort_trans;
    Float_t singleGausFWTM_normal;
    Float_t singleGausFWTM_trans;
    Float_t singleGaus_normal;
    Float_t singleGaus_trans;
    std::set<Int_t> maskedChannels;
    ClassDef(TResults,7);
};

#endif /* TRESULTS_HH_ */
