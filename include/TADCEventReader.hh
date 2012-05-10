/*
 * TADCEvent.hh
 *
 *  Created on: 01.08.2011
 *      Author: Felix Bachmair
 */

#ifndef TADCEVENT_H_
#define TADCEVENT_H_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCluster.hh"
#include "TPlaneProperties.hh"
//#include "TPlane.hh"
#include "TEvent.hh"
class TADCEventReader {
public:
	TADCEventReader(std::string fileName,UInt_t runNumber);
	virtual ~TADCEventReader();
	bool GetNextEvent();
	virtual bool LoadEvent(UInt_t EventNumber);
	Long64_t GetEntries();
	bool isOK();
//    bool getCMNEvent_flag() const;
	bool isValidTrack();
	UInt_t getAdcValue(UInt_t det,UInt_t ch);
	Float_t getSignalInSigma(UInt_t det,UInt_t ch);
	Float_t getSignal(UInt_t det,UInt_t ch);
    UInt_t getCurrent_event() const;
    UChar_t getDet_ADC(UInt_t i, UInt_t j) const;
    UChar_t getDet_Channels(UInt_t i , UInt_t j) const;
    UInt_t getDet_NChannels(UInt_t i) const;
    Float_t getDet_PedMean(UInt_t i, UInt_t j) const;
    Float_t getDet_PedWidth(UInt_t i, UInt_t j) const;
    UShort_t getDia_ADC(UInt_t i) const;
    UInt_t getEvent_number() const;
    TTree *getPedTree() const;
    UInt_t getRun_number() const;
    Float_t getStore_threshold() const;
    UInt_t getVerbosity() const;
    bool getZeroDivisorEvent_flag() const;
    TTree *getTree() const;
    static std::string getStringForDetector(int i);
    TFile* getFile() const;
    std::string getFilePath();
    Float_t getPedestalMean(UInt_t det, UInt_t ch);
    Float_t getPedestalSigma(UInt_t det, UInt_t ch);
    TCluster getCluster(UInt_t det,UInt_t cl);
    TCluster getCluster(UInt_t plane,TPlaneProperties::enumCoordinate cor, UInt_t cl);
    UInt_t getClusterSize(UInt_t det,UInt_t cl);
    UInt_t getNClusters(UInt_t det);
    bool isSaturated(UInt_t det,UInt_t ch);
    void checkADC();
    UInt_t getNDiamondClusters();
    bool isInFiducialCut();
    bool isDetMasked();
    TEvent* getEvent();
    void setVerbosity(UInt_t verbosity);
    bool useForAlignment(){return this->bUseForAlignment;};
    bool useForSiliconAlignment(){return this->bUseForSiliconAlignment;};
    bool useForAnalysis(){return this->bUseForAnalysis;};
    TH1F* getEtaIntegral(UInt_t det);
private:
	void SetBranchAddresses();
	bool SetTree(std::string fileName);//TTree *tree);
	void initialiseTree();
	int hasTree();
	TObject* getTreeName();
	void LoadEtaDistributions(UInt_t runNumber);
private:
	std::string fileName;
	UInt_t run_number;
	UInt_t event_number;
	Float_t store_threshold;
//	bool CMNEvent_flag;
	bool ZeroDivisorEvent_flag;
	UInt_t Det_NChannels[9];
	UChar_t Det_Channels[9][256];
	UChar_t Det_ADC[8][256];
	UShort_t Dia_ADC[128];
	Float_t Det_PedMean[9][256];
	Float_t Det_PedWidth[9][256];
	Float_t pedestalMean[9][256];
	Float_t pedestalSigma[9][256];
	//TCluster::vecvecTCluster* pVecvecCluster;
	TEvent *pEvent;
	bool bIsDetMasked;
	bool hasValidSiliconTrack;
	UInt_t nDiamondClusters;
	bool bIsInFiducialCut;
	//vector<bool> maskedDiaClusters;
	bool bUseForAlignment;
	bool bUseForAnalysis;
	bool bUseForSiliconAlignment;
private:
	TH1F *hEtaIntegral[9];
	bool bEtaIntegrals;
private:
//is that needed?
	TFile *file;
    TTree *tree;
    UInt_t current_event;
    TSystem* sys;

private:
	UInt_t verbosity;
};

#endif /* TADCEVENT_H_ */
