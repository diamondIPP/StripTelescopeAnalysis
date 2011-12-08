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
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"

#include "TADCEventReader.hh"
#include "TCluster.hh"
#include "TSettings.class.hh"

class TSelectionClass {
public:
	TSelectionClass(TSettings *settings);
	virtual ~TSelectionClass();
	void MakeSelection();
	void MakeSelection(UInt_t nEvents);
private:
	void setBranchAdressess();
	bool createSelectionTree(int nEvents);
	void resetVariables();
	void setVariables();
	bool checkDetMasked(UInt_t det);
	bool checkDetMasked(UInt_t det,UInt_t cl);
	TSettings *settings;
	TSystem *sys;
	TADCEventReader *eventReader;
	TFile* selectionFile;
	TTree* selectionTree;
	bool createdNewTree;
	bool createdNewFile;
	stringstream runString;
	stringstream rawfilepath;
	stringstream pedestalfilepath;
	stringstream clusterfilepath;
	HistogrammSaver *histSaver;
	UInt_t verbosity;
private:
	UInt_t nEvent;
	bool isDetMasked;//one of the Silicon Planes contains a Cluster with a masked channel
	vector<bool> isDiaMasked;//thediamond plane contains a cluster wit a masked channel (size of nDiamondHits)
	UInt_t nDiamondHits; //number of clusters in diamond plane;
	bool hasValidSiliconTrack; //One and only one cluster in each silicon plane;
	bool isInFiducialCut; //if hasValidSiliconTrack avarage of x and y of all planes is in fidcut region
	bool createdTree;
};

#endif /* TSELECTIONCLASS_HH_ */