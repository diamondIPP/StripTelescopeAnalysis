/*
 * THTMLCluster.hh
 *
 *  Created on: Apr 30, 2012
 *      Author: bachmair
 */

#ifndef THTMLCLUSTER_HH_
#define THTMLCLUSTER_HH_

#include "THTMLGenerator.hh"
#include "TADCEventReader.hh"
#include "TPlaneProperties.hh"
class THTMLCluster: public THTMLGenerator {
public:
	THTMLCluster(TSettings *settings);
	virtual ~THTMLCluster();
public:
	void createTableOfCuts();
};

#endif /* THTMLCLUSTER_HH_ */
