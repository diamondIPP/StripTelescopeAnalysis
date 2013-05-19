/**
 * @file TDiamondPattern.hh
 *
 * @date Jan 9, 2013
 * @author bachmair
 * @description
 */

#ifndef TDIAMONDPATTERN_HH_
#define TDIAMONDPATTERN_HH_

//C++ standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <utility> //pair
#include "TROOT.h"
#include "TPlaneProperties.hh"


/*
 *
 */
class TDiamondPattern:public TObject {
public:
	TDiamondPattern();
	virtual ~TDiamondPattern();
	void loadStandardPitchWidthSettings();
	bool addPattern(Float_t pitchWidth, Float_t startPosition, UInt_t firstChannel,UInt_t lastChannel);
	Float_t convertChannelToMetric(Float_t channel);
	Float_t convertMetricToChannel(Float_t metric);
	Float_t convertMetricToChannel(Float_t metric,UInt_t interval);
	Float_t getPitchWidth(UInt_t area);
	Int_t getNPatterns(){return nChannelsOfInterval.size();}
	void loadPitchWidthSettings(Float_t pitchWidth);
	void resetPattern();
	void clear(){resetPattern();}
	void Print();
	void showPatterns();
	UInt_t size(){return getNPatterns();}
	bool isStandardPitchWidth(){return bLoadedStandardPitchWidthSettings;}
	std::pair<Int_t,Int_t> getInterval(UInt_t pattern);
private:
	UInt_t getNIntervals(){return nChannelsOfInterval.size();}
	Float_t getChannelToMetric(UInt_t ch);
	void initialiseVector();
	std::vector<Float_t> channelToMetricConversion;
	std::vector<Float_t> beginOfInterval;
	std::vector<Float_t> endOfInterval;
	std::vector<Float_t> firstChannelOfInterval;
	std::vector<Float_t> nChannelsOfInterval;
	std::vector<Float_t> pitchWidth;
	bool bLoadedStandardPitchWidthSettings;
	Float_t standardPW;
    ClassDef(TDiamondPattern,1);
};

#endif /* TDIAMONDPATTERN_HH_ */
