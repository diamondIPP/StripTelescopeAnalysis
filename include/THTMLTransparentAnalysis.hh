/*
 * THTMLTransparentAnalysis.h
 *
 *  Created on: Jul 5, 2012
 *      Author: bachmair
 */

#ifndef THTMLTRANSPARENTANALYSIS_H_
#define THTMLTRANSPARENTANALYSIS_H_

#include "THTMLGenerator.hh"
#include "TSettings.class.hh"

class THTMLTransparentAnalysis: public THTMLGenerator {
public:
	
	THTMLTransparentAnalysis(TSettings *settings) ;
	virtual ~THTMLTransparentAnalysis();
	void createContent();
	void createPulseHeightPlots(vector<pair <UInt_t,Float_t> > meanPulseHeigths);
	void createResolutionPlots();
	void createEtaPlots();

private:
	UInt_t subjectDetector, subjectPlane;
};

#endif /* THTMLTRANSPARENTANALYSIS_H_ */
