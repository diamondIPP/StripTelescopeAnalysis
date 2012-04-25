/*
 * LandauGaussFit.hh
 *
 *  Created on: Apr 24, 2012
 *      Author: bachmair
 */

#ifndef LANDAUGAUSSFIT_HH_
#define LANDAUGAUSSFIT_HH_
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
class LandauGaussFit {
public:
	LandauGaussFit();
	virtual ~LandauGaussFit();
	TF1*  doLandauGaussFit(TH1F* histo,bool verbose=false);
private:
	TF1* langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues,
			Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams,
			Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF,bool verbose=false);
	Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);
//	Double_t langaufun(Double_t *x, Double_t *par);
};
#endif /* LANDAUGAUSSFIT_HH_ */

