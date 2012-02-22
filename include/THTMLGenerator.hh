/*
 * THTMLGenerator.hh
 *
 *  Created on: Feb 15, 2012
 *      Author: bachmair
 */

#ifndef THTMLGENERATOR_HH_
#define THTMLGENERATOR_HH_

#include <iostream>
#include <fstream>

#include "TSystem.h"
#include "TDatime.h"
#include "TSettings.class.hh"
#include "TDetectorAlignment.hh"

class THTMLGenerator {
public:
	THTMLGenerator(TSettings *settings);
	virtual ~THTMLGenerator();
	void generateHTMLFile();

private:
	void createAlignmentSummary(TDetectorAlignment *alignment);
	void generatorHTMLHeader();
	void generateHTMLTail();
	TSettings* settings;
	UInt_t verbosity;
	TSystem *sys;
	std::ofstream html_summary;
    TDatime dateandtime;
};

#endif /* THTMLGENERATOR_HH_ */
