/*
 * TRawEventSaver.cpp
 *
 *  Created on: 09.11.2011
 *      Author: bachmair
 */

#include "../include/TRawEventSaver.hh"

TRawEventSaver::TRawEventSaver(TSettings *settings){//unsigned int RunNumber, std::string RunDescription){
	// TODO Auto-generated constructor stub
	this->settings=settings;
	this->runNumber=settings->getRunNumber();
	rawEventReader = new TRawEventReader(settings);//this->runNumber);
	sys = gSystem;
	if(!settings->existsDirectory(settings->getRelativeOuputPath()))
	  sys->MakeDirectory(settings->getRelativeOuputPath().c_str());;
	sys->cd(settings->getRelativeOuputPath().c_str());

	rawfilepath.str("");
//	rawfilepath<<settings->getRelativeOuputPath().c_str();
	rawfilepath<<"rawData."<<runNumber<<".root";
	cout<<sys->pwd()<<endl;
	cout<<rawfilepath.str()<<endl;
	this->rawFile=TFile::Open(rawfilepath.str().c_str());
	cout<<"rawFile"<<rawFile<<"\t"<<flush;
	if(rawFile==0){
		createdNewFile=true;
		this->rawFile= new TFile(rawfilepath.str().c_str(),"recreate");
		cout<<"created new file"<<endl;
	}
	else
		createdNewFile=false;
	treeDescription<<"Raw Data of run "<<runNumber;
	rawFile->GetObject("rawTree",rawTree);
	if(rawTree==NULL){
		this->rawTree=new TTree("rawTree",treeDescription.str().c_str());
		createdNewTree=true;
	}
	else
		createdNewTree=false;
	sys->cd("..");
//
//
}

TRawEventSaver::~TRawEventSaver() {
	// TODO Auto-generated destructor stub
	cout<<"closing Files...,"<<endl;
	rawTree->Write("rawTree");
	rawFile->Close();
	sys->cd("..");
}


void TRawEventSaver::saveEvents(int nEvents){
	if(settings==0){
		cerr<<"SETTINGS Not initialized...";
		settings = new TSettings();
	}
	if (treeExists(nEvents)){
		cout<<"Tree has enough entries..."<<endl;
		return;
	}
	else{
		rawFile->Close();
		this->rawFile= new TFile(rawfilepath.str().c_str(),"recreate");
		this->rawTree=new TTree("rawTree",treeDescription.str().c_str());
		cout<<"Save Events to tree"<<endl;
		rawTree->Reset();
		this->setBranches();
		cout<<endl;
		for (int i=0;i<nEvents;i++){
			showStatusBar(i,nEvents,100);
			int  suceed = rawEventReader->ReadRawEvent(i,true);
			if(suceed<0){
			  cout<<"could not open file break"<<endl;
			  exit(-1);
			}
			loadEvent();
			eventNumber=i;
			//if(i%1000==0)		cout<<"Event:"<<i<<" \t\t"<<(int)Dia_ADC[1]<<endl;
			rawTree->Fill();
		}
	}
	cout<<"\nDONE."<<endl;
}

void TRawEventSaver::setBranches(){
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

void TRawEventSaver::loadEvent(){
	UInt_t diaInput=(UInt_t)settings->getDia_input();
	for(Int_t det=0;det<8;det++){
		for(Int_t ch=0; ch<256; ch++){
			Det_ADC[det][ch]=(UChar_t)rawEventReader->getPlane(det,diaInput).ADC_values[ch];
			//if(det==8)cout<<Det_ADC[det][ch]<<"\t"<<flush;
		}
	}
	for(UInt_t ch=0;ch<TPlaneProperties::getNChannelsDiamond();ch++){
		Dia_ADC[ch]=0;
	}
	for(UInt_t vaChNo=0;vaChNo<TPlaneProperties::getNChannelsDiamond();vaChNo++){
		//cout<<"Load diamond: "<<ch<<" -> "<<settings->getDetChannelNo(ch)<<endl;
		UInt_t detChNo= settings->getDetChannelNo(vaChNo);
		UShort_t value=(UShort_t)rawEventReader->getPlane(TPlaneProperties::getDetDiamond(),diaInput).ADC_values[vaChNo];
		if(detChNo<TPlaneProperties::getNChannelsDiamond())
			Dia_ADC[detChNo]=value;
		//if(ch==0)cout<<Dia_ADC[ch]<<" ";
	}
}

bool TRawEventSaver::treeExists(int nEvents){
	bool value=false;
	if(!this->createdNewFile&&!this->createdNewTree)
		if(this->rawTree->GetEntries()>=nEvents)
			value =true;
	cout<<"created new File = "<<createdNewFile<<" \tcreated new Tree: "<<createdNewTree<<" \tEntries:"<<rawTree->GetEntries()<<"of "<<nEvents<<"needed Entries, so treeExists: "<<value<<endl;
	return value;
}

void TRawEventSaver::showStatusBar(int nEvent,int nEvents,int updateIntervall,bool show){
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
	if(nEvent+1>=nEvents)cout<<endl;
}

