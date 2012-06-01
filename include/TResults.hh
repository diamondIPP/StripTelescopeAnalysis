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
#include "TSettings.class.hh"
#include <TDatime.h>

class TResults: public TNamed {
public:
  TResults(UInt_t runnumber=0);
  TResults(TSettings *settings);
  virtual ~TResults();
  void saveResults();
  void openResults(TSettings *Settings);
  TDatime getLastUpdateDate(){return lastUpdate;};
  void Print();

private:
  int runNumber;
//  TSettings Settings;
  TDatime lastUpdate;
  std::string path;
  UInt_t runnumber;
  ClassDef(TResults,1);
};

#endif /* TRESULTS_HH_ */
