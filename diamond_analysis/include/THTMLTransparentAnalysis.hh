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
	void createPulseHeightPlots(vector<vector <Float_t> > vecMeanPulseHeigths, vector<vector <Float_t> > vecMPPulseHeigths);
	void createResolutionPlots(vector<vector <pair <Float_t,Float_t> > > resolutions);
	void createEtaPlots();
	void createEtaIntegrals();

private:
	UInt_t subjectDetector, subjectPlane;
};

#endif /* THTMLTRANSPARENTANALYSIS_H_ */
