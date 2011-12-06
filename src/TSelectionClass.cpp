/*
 * TSelectionClass.cpp
 *
 *  Created on: 02.12.2011
 *      Author: bachmair
 */

#include "../include/TSelectionClass.hh"

TSelectionClass::TSelectionClass(TSettings* settings) {
	// TODO Auto-generated constructor stub
	cout<<"**********************************************************"<<endl;
	cout<<"************TSelectionClass::TSelectionClass**************"<<endl;
	cout<<"**********************************************************"<<endl;
	this->settings=settings;
	cout<<settings->getRunNumber()<<endl;

	// TODO Auto-generated constructor stub
	sys = gSystem;
	stringstream  runString;
	runString.str("");
	runString<<settings->getRunNumber();
	sys->MakeDirectory(runString.str().c_str());

	createdNewTree=false;
	createdNewFile=false;
	selectionTree=NULL;
	selectionFile=NULL;

	sys->cd(runString.str().c_str());
	rawfilepath<<"rawData."<<settings->getRunNumber()<<".root";
	pedestalfilepath.str("");
	pedestalfilepath<<"pedestalData."<<settings->getRunNumber()<<".root";
	clusterfilepath<<"clusterData."<<settings->getRunNumber()<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<clusterfilepath.str()<<endl;
	eventReader=new TADCEventReader(clusterfilepath.str());
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("selctions");
	sys->cd("selections");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
	sys->cd("..");
	verbosity=0;

}

TSelectionClass::~TSelectionClass() {
	// TODO Auto-generated destructor stub
	cout<<"\n\nClosing TSelectionClass"<<endl;
	delete eventReader;
	delete histSaver;
}

void TSelectionClass::MakeSelection()
{
	MakeSelection(eventReader->GetEntries());
}



void TSelectionClass::MakeSelection(UInt_t nEvents)
{
	cout<<"start selection with "<<nEvents<<" Events"<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		resetVariables();
		setVariables();

	}
}
bool TSelectionClass::createSelectionTree(int nEvents)
{
	bool createdNewFile=false;
	bool createdNewTree=false;
	stringstream selectionfilepath;
	selectionfilepath<<sys->pwd();
	selectionfilepath<<"/clusterData."<<settings->getRunNumber()<<".root";
	cout<<"Try to open \""<<selectionfilepath.str()<<"\""<<endl;
	selectionFile=new TFile(selectionfilepath.str().c_str(),"READ");
	if(selectionFile->IsZombie()){
		cout<<"clusterfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		selectionFile= new TFile(clusterfilepath.str().c_str(),"CREATE");
		selectionFile->cd();
	}
	else{
		createdNewFile=false;
		cout<<"File exists"<<endl;
	}
	selectionFile->cd();
	stringstream treeDescription;
	treeDescription<<"Selection Data of run "<<settings->getRunNumber();
	selectionFile->GetObject("selectionTree",selectionFile);
	if(selectionTree!=NULL){
		cout<<"File and Tree Exists... \t"<<flush;
		if(selectionTree->GetEntries()>=nEvents){
			createdNewTree=false;
			selectionTree->GetEvent(0);
				return false;
		}
		else{
			selectionTree->Delete();
			selectionTree=NULL;
		}
	}
	if(selectionTree==NULL){
		selectionFile->Close();
		selectionFile=new TFile(selectionfilepath.str().c_str(),"RECREATE");
		this->selectionTree=new TTree("selectionTree",treeDescription.str().c_str());
		createdNewTree=true;
		cout<<"\n\n\n***************************************************************\n";
		cout<<"***************************************************************\n";
		cout<<"there exists no tree:\'clusterTree\"\tcreate new one."<<selectionTree<<"\n";
		cout<<"***************************************************************\n";
		cout<<"***************************************************************\n"<<endl;
	}

	return createdNewTree;
}



void TSelectionClass::resetVariables(){

	isDetMasked = false;//one of the Silicon Planes contains a Cluster with a masked channel
	isDiaMasked.clear();
	nDiamondHits=0;
	hasValidSiliconTrack=true;; //One and only one cluster in each silicon plane;
}

void TSelectionClass::setVariables(){
	for(UInt_t det=0;det<8;det++){
//		cout<<(eventReader->getNClusters(det)==1);
		hasValidSiliconTrack=hasValidSiliconTrack&&(eventReader->getNClusters(det)==1);
		isDetMasked+=checkDetMasked(det);
	}
//	cout<<" "<<eventReader->getNClusters(8)<<":";
	for(UInt_t cl=0;cl<eventReader->getNClusters(8);cl++){
		isDiaMasked.push_back(checkDetMasked(8,cl));
//		cout<<isDiaMasked[cl];
	}
	cout<<endl;
//		bool isDetMasked;//one of the Silicon Planes contains a Cluster with a masked channel
//		vector<bool> isDiaMasked;//thediamond plane contains a cluster wit a masked channel (size of nDiamondHits)
//		UInt_t nDiamondHits; //number of clusters in diamond plane;
//		bool hasValidSiliconTrack; //One and only one cluster in each silicon plane;
//		bool isInFiducialCut; //if hasValidSiliconTrack avarage of x and y of all planes is in fidcut region
}



bool TSelectionClass::checkDetMasked(UInt_t det){
	bool isMasked=false;

	for(UInt_t cl=0;cl<eventReader->getNClusters(det);cl++){
		isMasked=isMasked||checkDetMasked(det,cl);
	}
	return isMasked;
}

bool TSelectionClass::checkDetMasked(UInt_t det,UInt_t cl){
	bool isMasked=false;

	if(cl<eventReader->getNClusters(det)){
		TCluster cluster = eventReader->getCluster(det,cl);
		for(UInt_t ch = cluster.getMinChannelNumber(); ch<=cluster.getMaxChannelNumber();ch++){
			isMasked=isMasked||settings->getDet_channel_screen(det).isScreened(ch);
		}
	}
	else
		return true;
	return isMasked;

	return false;
}
