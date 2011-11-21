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
class TADCEventReader {
public:
	TADCEventReader(std::string fileName);
	virtual ~TADCEventReader();
	bool GetNextEvent();
	bool GetEvent(UInt_t EventNumber);
	Long64_t GetEntries();
	bool isOK();
//    bool getCMNEvent_flag() const;
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
    static std::string getStringForPlane(int i);
    TFile* getFile() const;
    std::string getFilePath();
    Float_t getPedestalMean(UInt_t det, UInt_t ch);
    Float_t getPedestalSigma(UInt_t det, UInt_t ch);

private:
	void SetBranchAddresses();
	bool SetTree(std::string fileName);//TTree *tree);
	void initialiseTree();
	int hasTree();
	TObject* getTreeName();
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
private:
//is that needed?
	TFile *file;
    TTree *tree;
    UInt_t current_event;

private:
	UInt_t verbosity;
};

#endif /* TADCEVENT_H_ */
