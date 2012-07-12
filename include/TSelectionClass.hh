/*
 * TSelectionClass.hh
 *
 *  Created on: 02.12.2011
 *      Author: bachmair
 */

#ifndef TSELECTIONCLASS_HH_
#define TSELECTIONCLASS_HH_
#include "TSystem.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TPie.h"
#include "THTMLSelection.hh"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"

#include "TADCEventReader.hh"
#include "TCluster.hh"
#include "TSettings.class.hh"
#include "TResults.hh"

class TSelectionClass {
public:
	TSelectionClass(TSettings *settings);
	virtual ~TSelectionClass();
	void MakeSelection();
	void MakeSelection(UInt_t nEvents);
	void SetResults(TResults *results){this->results=results;};
private:
	TResults *results;
	void setBranchAdressess();
	bool createSelectionTree(int nEvents);
	void createCutFlowDiagramm();
	void resetVariables();
	void setVariables();
	bool isSaturated(UInt_t det,UInt_t cl=0);
	bool checkDetMasked(UInt_t det);
	bool checkDetMasked(UInt_t det,UInt_t cl);
	void initialiseHistos();
	void saveHistos();
	TSettings *settings;
	TSystem *sys;
	TADCEventReader *eventReader;
	TFile* selectionFile;
	TTree* selectionTree;
	bool createdNewTree;
	bool createdNewFile;
	bool isDiaSaturated;
	stringstream runString;
	stringstream rawfilepath;
	stringstream pedestalfilepath;
	stringstream clusterfilepath;
	HistogrammSaver *histSaver;
	UInt_t verbosity;
private:
	void findFiducialCut(TH2F* hFidCut);
	std::vector<Float_t >findFiducialCutIntervall(TH1D* hProj);
	THTMLSelection *htmlSelection;
	UInt_t nEvent;
	bool isDetMasked;//one of the Silicon Planes contains a Cluster with a masked channel
	deque<bool> isDiaMasked;//thediamond plane contains a cluster wit a masked channel (size of nDiamondHits)
	UInt_t nDiamondHits; //number of clusters in diamond plane;
	bool hasValidSiliconTrack; //One and only one cluster in each silicon plane;
	bool isInFiducialCut; //if hasValidSiliconTrack avarage of x and y of all planes is in fidcut region
	bool createdTree;
	bool useForAlignment;
	bool useForAnalysis;
	bool useForSiliconAlignment;
	Int_t nDiaClusterSize;
	UInt_t nToBigDiamondCluster;
	UInt_t nValidButMoreThanOneDiaCluster;
	UInt_t nValidSiliconNoDiamondHit;
	UInt_t nUseForAlignment;
	UInt_t nUseForSiliconAlignment;
	UInt_t nUseForAnalysis;
	UInt_t nNoValidSiliconTrack;
	UInt_t nValidSiliconAndDiamondCluster;
	UInt_t nValidDiamondTrack;
	UInt_t nValidSiliconTrack;
	UInt_t nSiliconTrackNotFiducialCut;
	UInt_t nEvents;

private:
	TH2F *hFiducialCutSilicon;
	TH2F *hFiducialCutSiliconDiamondHit;
};

#endif /* TSELECTIONCLASS_HH_ */
