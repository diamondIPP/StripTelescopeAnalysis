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
  virtual ~TResults();
  TResults(const TResults& rhs);
  void saveResults();
  void openResults(TSettings *Settings);
  TDatime getLastUpdateDate(){return lastUpdate;};
  void Print();
  void SetNoise(UInt_t det, Float_t detNoise);

private:
  void initialiseResults();
  void inheritOldResults(const TResults &rhs);
  int runNumber;
//  TSettings Settings;
  TDatime lastUpdate;
  std::string path;
  UInt_t runnumber;
//  TSettings *settings;
  std::vector<Float_t> seedSigma;
  std::vector<Float_t> hitSigma;
  std::vector<Float_t> noise;
  ClassDef(TResults,2);
};

#endif /* TRESULTS_HH_ */
