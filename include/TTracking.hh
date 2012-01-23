/*
 * TTracking.hh
 *
 *  Created on: Jan 23, 2012
 *      Author: bachmair
 */

#ifndef TTRACKING_HH_
#define TTRACKING_HH_

#include "TEventReader.hh"
#include "TTrack.hh"
#include "TDetectorAlignment.hh"
#include "TFile.h"

class TTracking: public TADCEventReader{
public:
	TTracking(std::string pathName, std::string alignmentName);
	virtual ~TTracking();
	bool LoadEvent(UInt_t eventNumber);
private:
	bool setAlignment(std::string alignmentName);
	TTrack *myTrack;
	TFile* alignmentFile;
	TDetectorAlignment* myAlignment;
};

#endif /* TTRACKING_HH_ */
