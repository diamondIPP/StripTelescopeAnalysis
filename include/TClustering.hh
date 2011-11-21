/*
 * TClustering.hh
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#ifndef TCLUSTERING_HH_
#define TCLUSTERING_HH_

using namespace std;

#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>

#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"

#include "TADCEventReader.hh"
class TClustering {
public:
	TClustering(int runNumber,int seedSigma,int hitSigma);
	virtual ~TClustering();
private:
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    int nEvent;
    int seedSigma;
    int hitSigma;
};

#endif /* TCLUSTERING_HH_ */
