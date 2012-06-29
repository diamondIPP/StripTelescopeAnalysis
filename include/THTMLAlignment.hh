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
	void setAlignment(TDetectorAlignment *alignment);
public:
	void createContent();

private:
	TDetectorAlignment* alignment;

};

#endif /* THTMLALIGNMENT_H_ */
