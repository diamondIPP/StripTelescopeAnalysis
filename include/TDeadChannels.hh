/*
 * TDeadChannels.hh
 *
 *  Created on: 18.11.2011
 *      Author: bachmair
 */

#ifndef TDEADCHANNELS_HH_
#define TDEADCHANNELS_HH_

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

#define N_DET_CHANNELS 256
using namespace std;

class TDeadChannels {
public:
	TDeadChannels(int runnumber,int seedSigma=10,int hitSigma=7);
	virtual ~TDeadChannels();
	void	doAnalysis(int nEvents=0);
private:
	void checkForDeadChannels();
	void initialiseHistos();
	void checkForSaturatedChannels();
	TH1F *hSaturatedChannels[8];
	TH1F *hSeedMap[8];
	TH1F *hClusterMap[8];
	TH1F* hNumberOfSeeds[8];
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
    int nEvent;
    int seedSigma;
    int hitSigma;
};

#endif /* TDEADCHANNELS_HH_ */
