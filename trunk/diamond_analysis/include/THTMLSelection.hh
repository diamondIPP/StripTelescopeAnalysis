/*
 * THTMLSelection.hh
 *
 *  Created on: May 16, 2012
 *      Author: bachmair
 */

#ifndef THTMLSELECTION_HH_
#define THTMLSELECTION_HH_

#include "THTMLGenerator.hh"

class THTMLSelection: public THTMLGenerator {
public:
	THTMLSelection(TSettings *settings);
	virtual ~THTMLSelection();
	void createFiducialCuts();
	void createCutFlowTable(std::vector<int> vecCutFlow);
	void createCutFlowGraph(std::string content);
};

#endif /* THTMLSELECTION_HH_ */
