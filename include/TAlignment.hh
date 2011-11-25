/*
 * TAlignment.hh
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#ifndef TALIGNMENT_HH_
#define TALIGNMENT_HH_

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
class TAlignment {
public:
	TAlignment(int runNumber);
	virtual ~TAlignment();
	void createVectors(UInt_t nEvents);
private:
	void initialiseHistos();
	void saveHistos();
	TADCEventReader* eventReader;
	HistogrammSaver* histSaver;
    TSystem* sys;
};

#endif /* TALIGNMENT_HH_ */
