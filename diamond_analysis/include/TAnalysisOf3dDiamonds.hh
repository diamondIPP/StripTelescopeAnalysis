/*
 * TAnalysisOf3dDiamonds.hh
 *
 *  Created on: Nov 20, 2012
 *      Author: bachmair, iain
 */

#ifndef TANALYSISOF3DDIAMONDS_HH_
#define TANALYSISOF3DDIAMONDS_HH_
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <deque>

#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "THTMLPedestal.hh"
#include "LandauGaussFit.hh"

#include "TADCEventReader.hh"
#include "TTracking.hh"
#include "TSettings.class.hh"
using namespace std;

class TAnalysisOf3dDiamonds {
public:
	TAnalysisOf3dDiamonds(TSettings *settings);
	virtual ~TAnalysisOf3dDiamonds();
	void	doAnalysis(UInt_t nEvents=0);
private:
	void initialiseHistos();
	void saveHistos();
	void analyseEvent();
private:
	TSettings *settings;
	HistogrammSaver *histSaver;
	TTracking* eventReader;
	UInt_t nEvent;
	vector<Float_t> vecXPredicted,vecYPredicted,vecXPredictedDiamondHit,vecYPredictedDiamondHit,
					vecPHDiamondHit,vecXPredictedDiamondHitFidCut,vecYPredictedDiamondHitFidCut,
					vecXFidCut,vecYFidCut,vecSpecialChannelMetricPos,vecChannel, vecChannelMetricPos;
	Int_t subjectPlane;
	vector<UInt_t> vecSilPlanes;

	TPositionPrediction *predictedPosition;
	TH1F* hChi2X;
	TH1F* hChi2Y;
	TH2F* hChi2XY;
	TH1F* hChi2;
	Int_t verbosity;

	TH2F* hValidSiliconAndOneDiamondHit;
	TH2F* hValidSiliconAndOneDiamondHitNotMasked;
	TH2F* hValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels;
	TH2F* hValidSiliconAndOneDiamondHitInOneArea;
	TH2F* hValidSiliconAndOneDiamondHitInSameAreaAndFidCut;
	TH2F* hAreaVsFidCut;

	UInt_t events;
	UInt_t noValidTrack;
	UInt_t tooHighChi2;
	UInt_t invalidPositionPrediction;
	UInt_t notOneDiaCluster;
	UInt_t notInFidCut;
	UInt_t saturatedCluster;
	UInt_t areaAndFidCutDoNotAgree;
	UInt_t hitNotInOneArea;
	UInt_t validEvents;
	UInt_t noDiaCluster;
	UInt_t maskedClusters;

};

#endif /* TANALYSISOF3DDIAMONDS_HH_ */
