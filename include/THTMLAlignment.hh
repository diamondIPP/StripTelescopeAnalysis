/*
 * THTMLAllignment.h
 *
 *  Created on: May 14, 2012
 *      Author: bachmair
 */

#ifndef THTMLALLIGNMENT_H_
#define THTMLALIGNMENT_H_

#include "THTMLGenerator.hh"

class THTMLAlignment: public THTMLGenerator {
public:
	THTMLAlignment(TSettings* settings);
	virtual ~THTMLAlignment();
	void setAlignment(TDetectorAlignment *alignment){this->alignment=alignment;};
public:
	void createContent();

private:
	void createOverviewTable();
	void createPostSiliconOverview();
	void createPostDiamondOverview();
	void createPreSiliconOverview();
	void createPreDiamondOverview();
	TDetectorAlignment* alignment;

};

#endif /* THTMLALIGNMENT_H_ */
