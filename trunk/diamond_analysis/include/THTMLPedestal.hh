/*
 * THTMLPedestal.hh
 *
 *  Created on: May 2, 2012
 *      Author: bachmair
 */

#ifndef THTMLPEDESTAL_HH_
#define THTMLPEDESTAL_HH_

#include "THTMLGenerator.hh"
#include "TADCEventReader.hh"
#include "TPlaneProperties.hh"

class THTMLPedestal: public THTMLGenerator {
public:
	THTMLPedestal(TSettings* settings);
	virtual ~THTMLPedestal();
public:
	void createPageContent();
	void createTableOfCuts();
	void createPedestalDistribution();
	void createBiggestHitMaps();
	void createNoiseDistribution();
	void createHitOrderSection();
	void createSaturatedChannels();
};

#endif /* THTMLPEDESTAL_HH_ */
