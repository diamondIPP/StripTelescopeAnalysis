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
#include <cstring>
#include <TObject.h>
#include <vector>
#include "TPlaneProperties.hh"
#include "TSettings.class.hh"
#include <TDatime.h>

class TResults: public TNamed {
public:
  TResults(UInt_t runnumber=0);
  TResults(TSettings *settings);
  TResults(const   TResults& rhs);//copy constructor
  TResults &operator=(const   TResults &src); //class assignment function
  virtual ~TResults();
  TResults(const TResults& rhs);
  void saveResults();
  void openResults(TSettings *Settings);
  void Print();

public:
  void setAlignment(TDetectorAlignment newAlignment);
  TDatime getLastUpdateDate(){return lastUpdate;};
  void SetNoise(UInt_t det, Float_t detNoise);
  UInt_t getAllEvents() const{return nAllEvents;}
  UInt_t getExactlyOneDiamondHit() const { return nExactlyOneDiamondHit;}
  UInt_t getMoreThanOneDiamondHit() const {return nMoreThanOneDiamondHit;}
  UInt_t getNoDiamondHit() const {return nNoDiamondHit;}
  UInt_t getNoSiliconHit() const {return nNoSiliconHit;}
  UInt_t getOneAndOnlyOneSiliconNotFiducialCut() const {return nOneAndOnlyOneSiliconNotFiducialCut;}
  UInt_t getValidSiliconTrack() const {return nValidSiliconTrack;}
  UInt_t getUseForAlignment() const {return nUseForAlignment;}
  UInt_t getUseForAnalysis() const {return nUseForAnalysis;}
  void setAllEvents(UInt_t allEvents){nAllEvents = allEvents;}
  void setExactlyOneDiamondHit(UInt_t exactlyOneDiamondHit){ nExactlyOneDiamondHit = exactlyOneDiamondHit;}
  void setMoreThanOneDiamondHit(UInt_t moreThanOneDiamondHit){ nMoreThanOneDiamondHit = moreThanOneDiamondHit;}
  void setNoDiamondHit(UInt_t noDiamondHit){ nNoDiamondHit = noDiamondHit;}
  void setNoSiliconHit(UInt_t noSiliconHit){ nNoSiliconHit = noSiliconHit;}
  void setOneAndOnlyOneSiliconNotFiducialCut(UInt_t oneAndOnlyOneSiliconNotFiducialCut){ nOneAndOnlyOneSiliconNotFiducialCut = oneAndOnlyOneSiliconNotFiducialCut;}
  void setValidSiliconTrack(UInt_t validSiliconTrack){ nValidSiliconTrack = validSiliconTrack;};
  void setUseForAlignment(UInt_t useForAlignment){ nUseForAlignment = useForAlignment;}
  void setUseForAnalysis(UInt_t useForAnalysis){nUseForAnalysis = useForAnalysis;};
public:
    void initialiseResults();
    void inheritOldResults(const TResults & rhs);
private:
    TDetectorAlignment alignment;
    int runNumber;
    //  TSettings Settings;
    TDatime lastUpdate;
    std::string path;
    UInt_t runnumber;
    //  TSettings *settings;
    std::vector<Float_t> seedSigma;
    std::vector<Float_t> hitSigma;
    std::vector<Float_t> noise;
    std::vector<Float_t> clusterSize;
    std::vector<Float_t> clusterSeedSize;
    std::vector<Float_t> nClusters;
    std::vector<Float_t> clusterPHmean;
    std::vector<Float_t> clusterPHwidth;
    UInt_t nAllEvents;
    UInt_t nNoSiliconHit;
    UInt_t nOneAndOnlyOneSiliconNotFiducialCut;
    UInt_t nValidSiliconTrack;
    UInt_t nNoDiamondHit;
    UInt_t nMoreThanOneDiamondHit;
    UInt_t nExactlyOneDiamondHit;
    UInt_t nUseForAlignment;
    UInt_t nUseForAnalysis;
    ClassDef(TResults,2);
};

#endif /* TRESULTS_HH_ */
