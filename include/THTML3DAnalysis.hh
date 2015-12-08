/*
 * THTML3DAnalysis.h
 *
 *  Created on: May 14, 2012
 *      Author: bachmair
 */

#ifndef THTML3DANALYSIS_H_
#define THTML3DANALYSIS_H_

#include "THTMLGenerator.hh"

class THTML3DAnalysis: public THTMLGenerator {
public:
	THTML3DAnalysis(TSettings* settings);
	virtual ~THTML3DAnalysis();
public:
	void createContent();

private:
	void createOverviewPlots();
	void createDistributionComparePlots();
	void createResolutionPlots();

};

#endif /* THTML3DANALYSIS_H_ */
