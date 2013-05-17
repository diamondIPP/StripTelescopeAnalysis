/**
 * @file TAnalysisOfAsymmetricEta.hh
 *
 * @date May 13, 2013
 * @author bachmair
 * @description
 */

#ifndef TANALYSISOFASYMMETRICETA_HH_
#define TANALYSISOFASYMMETRICETA_HH_

#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include "TSettings.class.hh"
#include "LandauGaussFit.hh"
#include "TProfile.h"
#include "TPolyMarker.h"
#include "TList.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"

#include "TSettings.class.hh"
#include "TCluster.hh"
#include "HistogrammSaver.class.hh"

/*
 *
 */
class TAnalysisOfAsymmetricEta {
public:
	TAnalysisOfAsymmetricEta(TSettings* settings);
	virtual ~TAnalysisOfAsymmetricEta();
	void setVectorOfCluster(vector<TCluster> vecClus);
	UInt_t analyse();
	Float_t getAlpha(){return alpha;}
	void setAlpha(Float_t newAlpha) {alpha=newAlpha;}
	void saveAsymmetricEtaPerArea(TH2F* histo, TString histName, Float_t alpha);
//	const vector<TCluster>& getClusters() const { return clusters; }
	void setClusters(vector<TCluster> clusters);// { this->clusters = clusters; }
	void Reset();
	UInt_t getDetector() const { return det; }
	void setDetector(UInt_t det) { this->det = det; }

private:
	UInt_t det;
	TH2F* hAsymmetricEta2D;
	void FillEtaDistribution(TH2F* hAsymmetricEta);
	TPolyMarker* FindPeaks(TH1F* histo, int nPeaks, Float_t sigma=4, TString option="Markov", Float_t threshold=.2);
	TPolyMarker* FindPeaks(UInt_t nTries,TH1F* histo, int nPeaks, Float_t sigma=4, TString option="Markov", Float_t threshold=.2);
	TH1F* getProjection();
	vector <TCluster> clusters;
	Float_t alpha;
	TSettings *settings;
	HistogrammSaver *histSaver;
	UInt_t maxTriesPeakfinding;
	UInt_t maxTriesAlpha;

	vector <Float_t> alphaValues;
	vector <Float_t> means;
	vector <Float_t> skewnesses;
	vector <Float_t> rightLefts;
	vector <Float_t> vecTries;
	int verbosity;
};

#endif /* TANALYSISOFASYMMETRICETA_HH_ */
