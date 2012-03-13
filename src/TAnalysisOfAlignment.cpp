/*
 * TAnalysisOfAlignment.cpp
 *
 *  Created on: Mar 2, 2012
 *      Author: bachmair
 */

#include "../include/TAnalysisOfAlignment.hh"

TAnalysisOfAlignment::TAnalysisOfAlignment(TSettings *settings) {
	// TODO Auto-generated constructor stub
	cout<<"\n\n\n\n**********************************************************"<<endl;
		cout<<"**********************************************************"<<endl;
		cout<<"*********TAnalysisOfClustering::TAnalysisOfAlignment*****"<<endl;
		cout<<"**********************************************************"<<endl;
		cout<<"**********************************************************\n\n\n"<<endl;
		// TODO Auto-generated constructor stub
		setSettings(settings);
		UInt_t runNumber=settings->getRunNumber();
		sys = gSystem;
		stringstream  runString;
		runString.str("");
		runString<<runNumber;
		sys->MakeDirectory(runString.str().c_str());

		sys->cd(runString.str().c_str());
		stringstream  filepath;
		filepath.str("");
		filepath<<"selectionData."<<runNumber<<".root";
		cout<<"currentPath: "<<sys->pwd()<<endl;
		cout<<filepath.str()<<endl;
		stringstream alignFile;
		alignFile<<"alignment."<<runNumber<<".root";
		eventReader=new TTracking(filepath.str(),alignFile.str(),runNumber);
		histSaver=new HistogrammSaver();
		sys->MakeDirectory("anaAlignmnet");
		sys->cd("anaAlignmnet");
		stringstream plotsPath;
		plotsPath<<sys->pwd()<<"/";
		histSaver->SetPlotsPath(plotsPath.str().c_str());
		histSaver->SetRunNumber(runNumber);
		sys->cd("..");
		initialiseHistos();
		cout<<"end initialise"<<endl;
		settings=0;
		verbosity=0;
}

TAnalysisOfAlignment::~TAnalysisOfAlignment() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TAnalysisOfAlignment::doAnalysis(UInt_t nEvents){
	if (nEvents==0||eventReader->GetEntries()<nEvents)
		nEvents=eventReader->GetEntries();
	cout<<"Do Analysis After Alignment...."<<endl;
	TH1F *histo = new TH1F("histo","histo",128,-0.501,0.501);
	for(nEvent=0;nEvent<nEvents;nEvent++){
		eventReader->LoadEvent(nEvent);
		if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
			continue;
		TCluster xClus = eventReader->getCluster(2,0);
		vector<UInt_t>refPlanes;
		refPlanes.push_back(0);
		refPlanes.push_back(2);
		refPlanes.push_back(3);
		TPositionPrediction* pred = eventReader->predictPosition(1,refPlanes,false);
		TH1F *hEtaIntegral=eventReader->getEtaIntegral(2);
		if(nEvent==0)histSaver->SaveHistogram(eventReader->getEtaIntegral(2));
		Float_t pos = eventReader->getStripXPositionOfCluster(1,xClus,pred->getPositionY(),TCluster::corEta,hEtaIntegral);
		Float_t stripX = eventReader->getStripXPosition(1,pos,TCluster::corEta);
		//cout<<"Event: "<<nEvent<<"\n\t "<<pred->getPositionX()<<"/"<<pred->getPositionY()<<"\n\t\t"<<xClus.getPosition(TCluster::corEta)<<" "<<pos<<" "<<stripX<<endl;
		UInt_t stripMiddle=(UInt_t) (stripX+0.5);
		Float_t delta = stripX-stripMiddle;
		histo->Fill(delta);
	}
	histSaver->SaveHistogram(histo);
	Int_t minBin = histo->GetMinimumBin();
	UInt_t nMinEntries =histo->GetBinContent(minBin);
	TH1F * histo2=new TH1F("histo2","histo2",128,-0.501,0.501);
	vector<UInt_t> vecEventNo;
	vecEventNo.clear();
	for(nEvent=0;nEvent<nEvents;nEvent++){
			eventReader->LoadEvent(nEvent);
			if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
				continue;
			TCluster xClus = eventReader->getCluster(2,0);
			vector<UInt_t>refPlanes;
			refPlanes.push_back(0);
			refPlanes.push_back(2);
			refPlanes.push_back(3);
			TPositionPrediction* pred = eventReader->predictPosition(1,refPlanes,false);
			TH1F *hEtaIntegral=eventReader->getEtaIntegral(2);
			if(nEvent==0)histSaver->SaveHistogram(eventReader->getEtaIntegral(2));
			Float_t pos = eventReader->getStripXPositionOfCluster(1,xClus,pred->getPositionY(),TCluster::corEta,hEtaIntegral);
			Float_t stripX = eventReader->getStripXPosition(1,pos,TCluster::corEta);
			UInt_t stripMiddle=(UInt_t) (stripX+0.5);
			Float_t delta = stripX-stripMiddle;
			Int_t bin =histo2->FindBin(delta);
			if(histo2->GetBinContent(bin)<nMinEntries){
				vecEventNo.push_back(nEvent);
				histo2->Fill(delta);
			}
	}
	histSaver->SaveHistogram(histo2);

	TH1F *hEta=new TH1F("hEtaDistribution","hEtaDistribution on det 2",128,0,1);
	for(UInt_t i=0;i<vecEventNo.size();i++){
		nEvent= vecEventNo.at(i);
		eventReader->LoadEvent(nEvent);

		if(!eventReader->useForAlignment()&&!eventReader->useForAnalysis())
			continue;
		Float_t eta = eventReader->getCluster(2,0).getEta();
		hEta->Fill(eta);

	}
	histSaver->SaveHistogram(hEta);
}
void TAnalysisOfAlignment::setSettings(TSettings *settings)
{
	this->settings=settings;
}

void TAnalysisOfAlignment::initialiseHistos()
{
}





