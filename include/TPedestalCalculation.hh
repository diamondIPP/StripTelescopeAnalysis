/*
 * PedestalCalculation.hh
 *
 *  Created on: 10.11.2011
 *      Author: bachmair
 */

#ifndef PEDESTALCALCULATION_HH_
#define PEDESTALCALCULATION_HH_

//C++ standard libraries
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <utility>
#include <deque>
#include <algorithm>
#include "TSystem.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TADCEventReader.hh"
#include "TRawEventSaver.hh"
#include "HistogrammSaver.class.hh"
#include <math.h>
#include "TMath.h"

using namespace std;

#define N_DET_CHANNELS 256
#define N_DIA_CHANNELS 128
class TPedestalCalculation {
public:
	TPedestalCalculation(TSettings *settings);
	//	TPedestalCalculation(int runNumber, int nEvents);
	virtual ~TPedestalCalculation();
	void calculateStartingPedestal(int nEvents);
	void calculateSlidingPedestals(UInt_t nEvents);
	static Float_t RoundFloat(Float_t value,UInt_t prec=2){if(value!=value)value=0;return (Float_t)((Long_t)(value*TMath::Power((Double_t)10,(Double_t)prec)+0.5))/TMath::Power((Double_t)10,(Double_t)prec);}
private:
	void calculateFirstPedestals(deque<Int_t> DetAdcQueue[8][N_DET_CHANNELS], deque<Int_t> DiaAdcQueue[N_DIA_CHANNELS],int maxSigma=7);
	pair <Float_t,Float_t> calculateFirstPedestalDet(int det,int ch, deque<Int_t> adcQueue, float mean, float sigma, int iterations=5,float maxSigma=7);
	pair <Float_t,Float_t> calculateFirstPedestalDia(int ch, deque<Int_t> adcQueue, float mean, float sigma, int iterations=5,float maxSigma=5);
	pair <Float_t,Float_t> calculateFirstPedestalDiaCMN(int ch, deque<Float_t> adcQueue, float mean, float sigma, int iterations=5,float maxSigma=5);
	pair <Float_t,Float_t> checkPedestalDet(int det, int ch,int maxSigma=7);
	pair <Float_t,Float_t> checkPedestalDia(int ch,int maxSigma=7);
	void printDiamond(UInt_t nChannel);
	bool createPedestalTree(int nEvents);
	void setBranchAdresses();
	void doCmNoiseCalculation();
	void fillFirstEventsAndMakeDiaDeque();
	void initialiseDeques();
	void updateDiamondPedestals();
	void updateSiliconPedestals();
	TADCEventReader* eventReader;
	TFile* pedestalFile;
	TTree* pedestalTree;
	bool createdNewTree;
	UInt_t nEvents;
	bool createdNewFile;
	bool doCMNCorrection;
	TSystem* sys;
	TSettings *settings;
	UInt_t runNumber;
	Float_t pedestalMean[9][N_DET_CHANNELS];
	Float_t  pedestalSigma[9][N_DET_CHANNELS];

	Float_t diaPedestalMean[N_DIA_CHANNELS];
	Float_t diaPedestalSigma[N_DIA_CHANNELS];

	Float_t diaPedestalMeanStartValues[N_DIA_CHANNELS];
	Float_t diaPedestalSigmaStartValues[N_DIA_CHANNELS];

	double sigmaValues[9][N_DET_CHANNELS];
	double meanValues[9][N_DET_CHANNELS];

	UInt_t slidingLength;
	deque<Int_t> detAdcValues[8][N_DET_CHANNELS];
	deque<Int_t> diaAdcValues[N_DIA_CHANNELS];
	deque<bool> detEventUsed[8][N_DET_CHANNELS];
	deque<bool> diaEventUsed[N_DIA_CHANNELS];

	UInt_t nEvent;
	Double_t detSUM[8][N_DET_CHANNELS];
	Double_t detSUM2[8][N_DET_CHANNELS];

	Double_t diaSUM[N_DIA_CHANNELS];
	Double_t diaSUM2[N_DIA_CHANNELS];
	Double_t diaSUMCmn[N_DIA_CHANNELS];
	Double_t diaSUM2Cmn[N_DIA_CHANNELS];


	int detEventsInSum[8][N_DET_CHANNELS];
	int diaEventsInSum[N_DIA_CHANNELS];
	int diaEventsInSumCMN[N_DIA_CHANNELS];


	Float_t diaPedestalMeanCMN[N_DIA_CHANNELS];
	Float_t diaPedestalSigmaCMN[N_DIA_CHANNELS];
	deque<Float_t> diaAdcValuesCMN[N_DIA_CHANNELS];
	deque<bool> diaEventUsedCMN[N_DIA_CHANNELS];
	Float_t cmNoise;

	//	stringstream rawfilepath;
	int MAXSDETSIGMA;
	int MAXDIASIGMA;
	HistogrammSaver *histSaver;
	TH1F* hCommonModeNoise;
	int printChannel;
	int verbosity;
	static Float_t getLowLimitPedestal(Float_t pedMean, Float_t pedSigma,Float_t maxSigma);
	static Float_t getHighLimitPedestal(Float_t pedMean, Float_t pedSigma,Float_t maxSigma);
};

#endif /* PEDESTALCALCULATION_HH_ */
