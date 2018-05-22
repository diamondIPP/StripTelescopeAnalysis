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
#include <map>

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"

using namespace std;
UInt_t eventNumber;
UInt_t runNumber;
Float_t chargeShareDia;
Float_t chargeShareSil;
UChar_t Det_ADC[8][256];
UShort_t Dia_ADC[128];
int verbosity = 5;
map<UInt_t,UInt_t> errorCounter;

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
//    rawTree->Branch("RunNumber",&runNumber,"RunNumber/i");
    rawTree->Branch("EventNumber",&eventNumber,"EventNumber/i");
}

void setBranchAdresses(TTree* tree){
    /*if(tree->FindBranch("RunNumber")){
        tree->SetBranchAddress("RunNumber",&runNumber);
        if(verbosity>3)cout<<"Set Branch \"RunNumber\""<<endl;
    }
    else if(tree->FindBranch("runNumber")){
        tree->SetBranchAddress("runNumber",&runNumber);
        if(verbosity>3)cout<<"Set Branch \"runNumber\""<<endl;
    }*/
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
    if(tree->FindBranch("DiaADC")){
        tree->SetBranchAddress("DiaADC",&Dia_ADC);
        if(verbosity>3)cout<<"Set Branch \"DiaADC\""<<endl;
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
void updateSilicon(){
    UShort_t max = 256;
    for (UInt_t det = 0; det < 8; det++){
        UInt_t startChannel = 0;
        UInt_t endChannel = 127;
        if(det == 2 || det == 6){
            startChannel = 127;
            endChannel = 0;
        }
        bool finished = false;
        Float_t alpha = chargeShareSil;
        UShort_t adc = 0;
        for(UInt_t ch=startChannel;!finished;){
            UShort_t measured_adc = Det_ADC[det][ch];
            Float_t real_adc = ((Float_t)measured_adc-(Float_t)adc*alpha)/(1-alpha);
            UShort_t newADC;
            if(real_adc> max -1)
                newADC = max -1;
            else
                newADC  = UShort_t(real_adc+.5);
            adc = newADC;
            Det_ADC[det][ch] = newADC;
            finished = (ch==endChannel);
            if(det==2||det == 6 )
                ch--;
            else
                ch++;
        }
    }
}

void updateDiamond(){
    Float_t alpha = chargeShareDia;
    UShort_t adc =  0;
    UShort_t max = 4096;
    TString out = "";
    out += TString::Format("%7d_%5.3f: ",eventNumber,alpha*100.);
    bool bPrint = false;
    for(UInt_t ch=0;ch<128;ch++){
        UShort_t measured_adc = Dia_ADC[ch];
        Float_t real_adc = ((Float_t)measured_adc-(Float_t)adc*alpha)/(1-alpha);
        UShort_t newADC;
        if(real_adc >= max-1)
            newADC = max-1;
        else
            newADC  = UShort_t(real_adc+.5);
        if (13<ch && ch< 30){
            out += TString::Format("[%2d, (%4d) %4d  --> %4d]  ",ch,adc,Dia_ADC[ch],newADC);
            if (TMath::Abs(Dia_ADC[ch]-newADC) > 1)
                bPrint = true;
        }
        adc = newADC;
        Dia_ADC[ch] = newADC;
    }
    if (bPrint&&false)
        cout<<out<<"\n"<<endl;
}

void LoopOverTree(TTree* inputTree,TTree* outputTree){
    setBranchAdresses(inputTree);
    UInt_t nEntries = inputTree->GetEntries();
    cout<<"Silicon Charge Share: "<<chargeShareSil<<"\tDiamond Charge Share: "<<chargeShareDia<<endl;
    //	char t; cin>>t;
    for (int i=0;i<inputTree->GetEntries();i++){

        showStatusBar(i,nEntries,100,0,0);
        inputTree->GetEntry(i);
        eventNumber =i;
        updateSilicon();
        updateDiamond();
        outputTree->Fill();
    }
}

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

int createAsymmetricEtaSample(int runNo=0,Float_t silCor=999,float diaCor=999,int corRunNo= -1){
    bool enter = (runNo==0);
    if (runNo==0){
        cout<<"Enter RunNumber: "<<flush;
        cin>>runNo;
    }
    runNumber = runNo;
    cout<<"The entered runnumber is "<<runNumber<<endl;

    if (silCor==999){
    	cout<<"\nHow much charge sharing do you want to activate in the Silicon (in %): "<<flush;
	    cin>>silCor;
    }
    chargeShareSil = silCor;
    chargeShareSil /=100;
    cout<<"There is charge sharing of "<<chargeShareSil<<" %%"<<endl;
    //for (UInt_t det = 0; det < chargeShareSil.size();i++){
    //    cout<<<<" * "<<det<<": "<<chargeShareSil[det]<<" %%"<<endl;
    //    chargeShareSil[det] /=100;
    //}
    if (diaCor==999){
	    cout<<"\nHow much charge sharing do you want to activate in the Diamond (in %): "<<flush;
	    cin>>diaCor;
    }
    chargeShareDia = diaCor;
	chargeShareDia /=100;
	cout<<"There is charge sharing of "<<chargeShareDia*100<<"%."<<endl;

	TString fileName = TString::Format("rawData.%05d.root",runNumber);
	cout<<"Reading file '"<<fileName<<"'"<<endl;

	TFile * file = (TFile*)TFile::Open(fileName);
	if (!file){
		cerr<<"File does not exists... EXIT"<<endl;
		return -1;
	}
	TTree* tree;
	file->GetObject("rawTree",tree);

	if(tree==NULL) {
		tree=(TTree*)getTreeName(file);
	}
	if (!tree){
		cerr<<"Tree does not exists... EXIT"<<endl;
		return -1;
	}
    cout<<"ChargeShareSil: "<<chargeShareSil<<"---> "<<(int)(chargeShareSil*1e4)<<endl;
    cout<<"ChargeShaReDia: "<<chargeShareDia<<"---> "<<(int)(chargeShareDia*1e4)<<endl;

	TString outputfileName = TString::Format("rawData.%05d-%05d-%05d.root",runNumber,(int)(chargeShareSil*1e4),(int)(chargeShareDia*1e4));
	cout<<outputfileName;
	if (!is_file_exist(outputfileName)){
        TFile* outputFile = new TFile(outputfileName,"RECREATE");
        outputFile->cd();
        TTree* outputTree =  new TTree("pedestalTree","pedestalTree");
        setBranches(outputTree);
        LoopOverTree(tree,outputTree);
        outputTree->Write("rawTree");
        outputFile->Close();
	}
	else{
	    cout<<"File already exists! "<<outputfileName<<endl;
	}
	TString command = ".! ln -s ";
	command += outputfileName;
    if (enter){
        cout<<"Enter Cor Runno: "<<flush;
        cin>>corRunNo;
    }
	if (corRunNo == -1)
	    corRunNo = runNumber*10+0;
    cout<<"Corrected Runno is "<<corRunNo<<endl;
	command += TString::Format(" rawData.%d.root",corRunNo);
	cout<<"Command: \""<<command<<"\""<<endl;
    TString linkFileName = TString::Format("rawData.%d.root",corRunNo);
    if (!is_file_exist(linkFileName))
        gROOT->ProcessLine(command);
    else
        cout<<"Link already exists: "<<linkFileName<<endl;
	return 1;
}



