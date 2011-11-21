/*
 * TEventReader.hh
 *
 *  Created on: 01.08.2011
 *      Author: Felix Bachmair
 */

#ifndef TEVENT_H_
#define TEVENT_H_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
class TEventReader {
public:
	TEventReader(std::string pathName, int runNumber);
	virtual ~TEventReader();
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

private:
	void SetBranchAddresses();
	bool SetTree(std::string fileName);//TTree *tree);
	void initialiseTree();
private:
	UInt_t run_number;
	UInt_t event_number;
	Float_t store_threshold;
//	bool CMNEvent_flag;
	bool ZeroDivisorEvent_flag;
	UChar_t Det_ADC[8][256];
	UShort_t Dia_ADC[256];
private:
//is that needed?
	TFile *file;
    TTree *tree;
    UInt_t current_event;

private:
	UInt_t verbosity;
};

#endif /* TEVENT_H_ */
