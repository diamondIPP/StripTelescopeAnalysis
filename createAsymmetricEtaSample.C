#include "TTree.h"
//C++ standard libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstring>
#include <deque>

#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"

using namespace std;
UInt_t eventNumber;
UInt_t runNumber;
Float_t chargeShare;
Float_t diaPedestalMeanCMN[128];
Float_t diaPedestalSigmaCMN[128];
Float_t diaPedestalMean[128];
Float_t diaPedestalMeanOld[128];
Float_t diaPedestalSigma[128];
UChar_t Det_ADC[8][256];
UShort_t Dia_ADC[128];
int verbosity = 5;

TObject* getTreeName(TFile* file){
	if(file==NULL) exit(-1);
	return 0;
}


void setBranches(TTree* rawTree){
	//rawTree->Branch("DetADC",&Det_ADC,"DetADC[8][256]/b");
	rawTree->Branch("D0X_ADC",&Det_ADC[0],"D0X_ADC[256]/b");
	rawTree->Branch("D0Y_ADC",&Det_ADC[1],"D0Y_ADC[256]/b");
	rawTree->Branch("D1X_ADC",&Det_ADC[2],"D1X_ADC[256]/b");
	rawTree->Branch("D1Y_ADC",&Det_ADC[3],"D1Y_ADC[256]/b");
	rawTree->Branch("D2X_ADC",&Det_ADC[4],"D2X_ADC[256]/b");
	rawTree->Branch("D2Y_ADC",&Det_ADC[5],"D2Y_ADC[256]/b");
	rawTree->Branch("D3X_ADC",&Det_ADC[6],"D3X_ADC[256]/b");
	rawTree->Branch("D3Y_ADC",&Det_ADC[7],"D3Y_ADC[256]/b");
	rawTree->Branch("DiaADC",&Dia_ADC,"DiaADC[128]/s");
	rawTree->Branch("RunNumber",&runNumber,"RunNumber/i");
	rawTree->Branch("EventNumber",&eventNumber,"EventNumber/i");
}

void setBranchAdresses(TTree* tree){
	if(tree->FindBranch("RunNumber")){
		tree->SetBranchAddress("RunNumber",&runNumber);
		if(verbosity>3)cout<<"Set Branch \"RunNumber\""<<endl;
	}
	else if(tree->FindBranch("runNumber")){
		tree->SetBranchAddress("runNumber",&runNumber);
		if(verbosity>3)cout<<"Set Branch \"runNumber\""<<endl;
	}
	if(tree->FindBranch("EventNumber")){
		tree->SetBranchAddress("EventNumber",&eventNumber);
		if(verbosity>3)cout<<"Set Branch \"EventNumber\""<<endl;
	}
	if(tree->FindBranch("D0X_ADC")){
			tree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
			if(verbosity>3)cout<<"Set Branch \"D0X_ADC\""<<endl;
		}
		if(tree->FindBranch("D0Y_ADC")){
			tree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
			if(verbosity>3)cout<<"Set Branch \"D0Y_ADC\""<<endl;
		}
		if(tree->FindBranch("D1X_ADC")){
			tree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
			if(verbosity>3)cout<<"Set Branch \"D1X_ADC\""<<endl;
		}
		if(tree->FindBranch("D1Y_ADC")){
			tree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
			if(verbosity>3)cout<<"Set Branch \"D1Y_ADC\""<<endl;
		}
		if(tree->FindBranch("D2X_ADC")){
			tree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
			if(verbosity>3)cout<<"Set Branch \"D2X_ADC\""<<endl;
		}
		if(tree->FindBranch("D2Y_ADC")){
			tree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
			if(verbosity>3)cout<<"Set Branch \"D2Y_ADC\""<<endl;
		}
		if(tree->FindBranch("D3X_ADC")){
			tree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
			if(verbosity>3)cout<<"Set Branch \"D3X_ADC\""<<endl;
		}
		if(tree->FindBranch("D3Y_ADC")){
			tree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
			if(verbosity>3)cout<<"Set Branch \"D3Y_ADC\""<<endl;
		}
		//tree->SetBranchAddress("Dia_ADC",&Dia_ADC);
		if(tree->FindBranch("DiaADC")){
			tree->SetBranchAddress("DiaADC",&Dia_ADC);
			if(verbosity>3)cout<<"Set Branch \"DiaADC\""<<endl;
		}

//		if(tree->FindBranch("Dia_PedMean")){
//			tree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
//			if(verbosity>3)cout<<"Set Branch \"Dia_PedMean\""<<endl;
//		}
//		if(tree->FindBranch("Dia_PedWidth")){
//			tree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
//			if(verbosity>3)cout<<"Set Branch \"Dia_PedWidth\""<<endl;
//		}
//		if(tree->FindBranch("PedestalMean")){
//			tree->SetBranchAddress("PedestalMean",&pedestalMean);
//			if(verbosity>3)cout<<"Set Branch \"PedestalMean\""<<endl;
//		}
//		if(tree->FindBranch("PedestalSigma")){
//			tree->SetBranchAddress("PedestalSigma",&pedestalSigma);
//			if(verbosity>3)cout<<"Set Branch \"PedestalSigma\""<<endl;
//		}
		if(tree->FindBranch("diaPedestalMean")){
			tree->SetBranchAddress("diaPedestalMean",&diaPedestalMean);
			if(verbosity>3)cout<<"Set Branch \"diaPedestalMean\""<<endl;
		}
		if(tree->FindBranch("diaPedestalSigma")){
			tree->SetBranchAddress("diaPedestalSigma",&diaPedestalSigma);
			if(verbosity>3)cout<<"Set Branch \"diaPedestalSigma\""<<endl;
		}
		if(tree->FindBranch("diaPedestalMeanCMN")){
			tree->SetBranchAddress("diaPedestalMeanCMN",&diaPedestalMeanCMN);
			if(verbosity>3)cout<<"Set Branch \"diaPedestalMeanCMN\""<<endl;
		}
		if(tree->FindBranch("diaPedestalSigmaCMN")){
			tree->SetBranchAddress("diaPedestalSigmaCMN",&diaPedestalSigmaCMN);
			if(verbosity>3)cout<<"Set Branch \"diaPedestalSigmaCMN\""<<endl;
		}
}
void showStatusBar(int nEvent,int nEvents,int updateIntervall,bool show,bool makeNewLine){
	if(nEvent+1>=nEvents)nEvent++;
	cout.precision(3);
	int percentageLength=50;
	if(nEvent%(int)updateIntervall==0||nEvent>=nEvents-1||show){
		double percentage = (double)(nEvent)/(double)nEvents*(double)100;
		cout<<"\rfinished with "<<setw(8)<<nEvent<<" of "<<setw(10)<<nEvents<<": "<<setw(6)<<std::setprecision(2)<<fixed<<percentage<<"%\t\tSTATUS:\t\t";
		for(int i=0;i<percentageLength;i++)
			if (i*10<percentage*(double)percentageLength/(double)10)cout<<"%";
			else cout<<"_";
		cout<<" "<<flush;
	}
	if(makeNewLine&&nEvent+1>=nEvents)cout<<endl;
}


void updateDiamond(){
	Float_t share =0;
	for(UInt_t ch=0;ch<128;ch++){
		UShort_t adc = Dia_ADC[ch];
		Float_t pedMean = diaPedestalMean[ch];
        if (pedMean<1000){
            cout<<"\ntaking old pedValue for: \t";
			cout<<eventNumber<<" "<<ch<<" "<<adc<<", "<<pedMean<<" = "<<", "<< share<<" --> "<<endl;
            pedMean = diaPedestalMeanOld[ch];
        }
        else diaPedestalMeanOld[ch] = pedMean;
		Float_t signal = (Float_t)adc-pedMean;
		//if(share>10)
        if(eventNumber >15637 &&eventNumber<15640&&false)
			cout<<eventNumber<<" "<<ch<<" "<<adc<<", "<<pedMean<<" = "<<signal<<", "<< share<<" --> ";
		adc += share;
		if(adc>4095)
            adc = 4095;
		if(eventNumber>15637 && eventNumber < 15640&&false)
				cout<<adc<<endl;
		Dia_ADC[ch] = adc;
		signal = (Float_t)adc-pedMean;
		share = signal * chargeShare;
	}
    if(eventNumber>15637 && eventNumber < 15640&&false){
        char t;
        cin>>t;
    }

}

void LoopOverTree(TTree* inputTree,TTree* outputTree){
    for(int i=0;i<128;i++)
        diaPedestalMeanOld[i]=0;
	setBranchAdresses(inputTree);
	UInt_t nEntries = inputTree->GetEntries();
	for (int i=0;i<inputTree->GetEntries();i++){
		showStatusBar(i,nEntries,100,0,0);
		inputTree->GetEntry(i);
		eventNumber =i;
		updateDiamond();
		outputTree->Fill();
	}
}


int createAsymmetricEtaSample(){
	cout<<"Enter RunNumber: "<<flush;
	cin>>runNumber;
	cout<<"The entered runnumber is "<<runNumber<<endl;

	cout<<"\nHow much charge sharing do you want to activate (in %): "<<flush;
	cin>>chargeShare;
	chargeShare /=100;
	cout<<"There is charge sharing of "<<chargeShare*100<<"%."<<endl;


	TString fileName = TString::Format("pedestalData.%05d.root",runNumber);
	cout<<"Reading file '"<<fileName<<"'"<<endl;

	TFile * file = (TFile*)TFile::Open(fileName);
	if (!file){
		cerr<<"File does not exists... EXIT"<<endl;
		return -1;
	}
	TTree* tree;
	file->GetObject("pedestalTree",tree);

	if(tree==NULL) {
		tree=(TTree*)getTreeName(file);
	}
	if (!tree){
		cerr<<"Tree does not exists... EXIT"<<endl;
		return -1;
	}
	TFile* outputFile = new TFile(TString::Format("rawData.%05d-%05d.root",runNumber,(int)(chargeShare*1e4)),"RECREATE");
	outputFile->cd();
	TTree* outputTree =  new TTree("pedestalTree","pedestalTree");
	setBranches(outputTree);
	LoopOverTree(tree,outputTree);
	outputTree->Write("rawTree");
	outputFile->Close();
	return 1;


}



