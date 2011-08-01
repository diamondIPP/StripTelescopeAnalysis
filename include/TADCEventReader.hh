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
using namespace std;
class TADCEventReader {
public:
	TADCEventReader(TTree* tree);
	virtual ~TADCEventReader();
	bool GetNextEvent();
	bool GetEvent(UInt_t EventNumber);
	Long64_t GetEntries();
private:
	void SetBranchAddresses();
	bool SetTree(TTree *tree);
	void initialiseTree();
public:
	UInt_t run_number;
	UInt_t event_number;
	Float_t store_threshold;
	bool CMNEvent_flag;
	bool ZeroDivisorEvent_flag;
	UInt_t Det_NChannels[9];
	UChar_t Det_Channels[9][256];
	UChar_t Det_ADC[8][256];
	UShort_t Dia_ADC[256];
	Float_t Det_PedMean[9][256];
	Float_t Det_PedWidth[9][256];
private:
//is that needed?
    TTree *PedTree;
    UInt_t current_event;

private:
	UInt_t verbosity;
};

#endif /* TADCEVENT_H_ */
