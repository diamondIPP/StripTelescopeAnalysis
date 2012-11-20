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
#include "TSettings.class.hh"
using namespace std;

class TRawEventSaver {
public:
	TRawEventSaver(TSettings *settings);//unsigned int RunNumber, std::string RunDescription = "");
	virtual ~TRawEventSaver();
	void setSettings(TSettings* set){settings=set;}
	void saveEvents(int nEvents);
	static void showStatusBar(int nEvent,int nEvents,int updateIntervall=100,bool show=false);
private:
	UInt_t runNumber;
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
    UInt_t eventNumber;
    TSettings *settings;
};

#endif /* TRAWEVENTSAVER_HH_ */
