/**
 * @file THTMLSelectionAnalysis.hh
 *
 * @date Feb 12, 2013
 * @author bachmair
 * @description
 */

#ifndef THTMLSELECTIONANALYSIS_HH_
#define THTMLSELECTIONANALYSIS_HH_

#include "THTMLGenerator.hh"
#include "TSettings.class.hh"

/*
 *
 */
class THTMLSelectionAnalysis: public THTMLGenerator {
public:
	THTMLSelectionAnalysis(TSettings *settings);
	virtual ~THTMLSelectionAnalysis();
	void addSelectionPlots();
	void addAreaPlots();
	void addFiducialCutPlots();
};

#endif /* THTMLSELECTIONANALYSIS_HH_ */
