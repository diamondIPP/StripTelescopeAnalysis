/*
 * TRawEventSaver.hh
 *
 *  Created on: 09.11.2011
 *      Author: bachmair
 */

#ifndef TRAWEVENTSAVER_HH_
#define TRAWEVENTSAVER_HH_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>

#include "TRawEventReader.hh"
#include "TFile.h"
#include "TTree.h"
using namespace std;

class TRawEventSaver {
public:
	TRawEventSaver(unsigned int RunNumber, std::string RunDescription = "");
	virtual ~TRawEventSaver();
	void saveEvents(int nEvents);
	static void showStatusBar(int nEvent,int nEvents,int updateIntervall=1000,bool show=false);
private:
	int runNumber;
	string runDesciption;
	TRawEventReader* rawEventReader;
	void setBranches();
	void loadEvent();
	bool treeExists(int nEvents);
	TFile *rawFile;
	TTree *rawTree;
    TSystem* sys;
	stringstream rawfilepath;
	stringstream treeDescription;
private:
    bool createdNewFile;
    bool createdNewTree;
    bool needToReloadEvents;
    UChar_t Det_ADC[8][256];
    UShort_t Dia_ADC[128];
    int eventNumber;
};

#endif /* TRAWEVENTSAVER_HH_ */
