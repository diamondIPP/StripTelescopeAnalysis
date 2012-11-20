/*
 * THTMLLandaus.hh
 *
 *  Created on: May 20, 2012
 *      Author: bachmair
 */

#ifndef THTMLLANDAUS_HH_
#define THTMLLANDAUS_HH_

#include "THTMLGenerator.hh"

class THTMLLandaus: public THTMLGenerator {
public:
	THTMLLandaus(TSettings *settings);
	virtual ~THTMLLandaus();
	void addLandauDiamondTable(vector<Float_t> vecHistoMeans,vector<Float_t> vecHistoMaxs,vector<Float_t> vecHistoGaus,vector<Float_t> vecHistoLandau);
	void addLandauDiamond(Float_t width,Float_t MP,Float_t area, Float_t GSigma);
};

#endif /* THTMLLANDAUS_HH_ */
