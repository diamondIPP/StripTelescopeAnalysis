/*
 * TSelectionClass.cpp
 *
 *  Created on: 02.12.2011
 *      Author: bachmair
 */

#include "../include/TSelectionClass.hh"

TSelectionClass::TSelectionClass(TSettings* settings) {
	// TODO Auto-generated constructor stub
	cout<<"\n\n\n	**********************************************************"<<endl;
	cout<<"************TSelectionClass::TSelectionClass**************"<<endl;
	cout<<"**********************************************************"<<endl;
	this->settings=settings;
	cout<<settings->getRunNumber()<<endl;

	// TODO Auto-generated constructor stub
	sys = gSystem;
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
	sys->MakeDirectory("selections");
	sys->cd("selections");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(settings->getRunNumber());
	sys->cd("..");
	cout<<"HISTSAVER:"<<sys->pwd()<<endl;

			verbosity=0;

	createdTree=false;

}

TSelectionClass::~TSelectionClass() {
	// TODO Auto-generated destructor stub
	cout<<"\n\nClosing TSelectionClass"<<endl;
	selectionFile->cd();
	if(selectionTree!=NULL&&this->createdTree){
		cout<<"CLOSING TREE"<<endl;
		cout<<eventReader->getTree()->GetName()<<" "<<clusterfilepath.str().c_str()<<endl;
		selectionTree->AddFriend(eventReader->getTree()->GetName(),clusterfilepath.str().c_str());
		cout<<"pedestalTree"<<" "<<pedestalfilepath.str().c_str()<<endl;
		selectionTree->AddFriend("pedestalTree",pedestalfilepath.str().c_str());
		cout<<"rawTree"<<" "<<rawfilepath.str().c_str()<<endl;
		selectionTree->AddFriend("rawTree",rawfilepath.str().c_str());
		cout<<"save selectionTree: "<<selectionTree->GetListOfFriends()->GetEntries()<<endl;
		selectionTree->Write();
	}
	selectionFile->Close();
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TSelectionClass::MakeSelection()
{
	MakeSelection(eventReader->GetEntries());
}



void TSelectionClass::MakeSelection(UInt_t nEvents)
{
	createdTree=createSelectionTree(nEvents);
	if(!createdTree) return;
	this->setBranchAdressess();
	cout<<"start selection with "<<nEvents<<" Events"<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		resetVariables();
		setVariables();
		selectionTree->Fill();
	}
}
bool TSelectionClass::createSelectionTree(int nEvents)
{
	cout<<"TSelectionClass::checkTree"<<endl;
	bool createdNewFile=false;
	bool createdNewTree=false;
	stringstream selectionfilepath;
	sys->cd(	runString.str().c_str());
	selectionfilepath<<sys->pwd();
	selectionfilepath<<"/selectionData."<<settings->getRunNumber()<<".root";
	cout<<"Try to open \""<<selectionfilepath.str()<<"\""<<endl;
	selectionFile=new TFile(selectionfilepath.str().c_str(),"READ");
	if(selectionFile->IsZombie()){
		cout<<"selectionfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		selectionFile= new TFile(selectionfilepath.str().c_str(),"CREATE");
		cout<<"DONE"<<flush;
		selectionFile->cd();
	}
	else{
		createdNewFile=false;
		cout<<"File exists"<<endl;
	}
	selectionFile->cd();
	cout<<"get Tree"<<endl;
	stringstream treeDescription;
	treeDescription<<"Selection Data of run "<<settings->getRunNumber();
	cout<<"get Tree2"<<endl;
	selectionFile->GetObject("selectionTree",selectionTree);
	cout<<"check Selection Tree:"<<selectionTree<<endl;
	if(selectionTree!=NULL){
		cout<<"File and Tree Exists... \t"<<flush;
		if(selectionTree->GetEntries()>=nEvents){
			createdNewTree=false;
			selectionTree->GetEvent(0);
				return false;
		}
		else{
			cout<<"selectionTree.events !- nEvents"<<flush;
			selectionTree->Delete();
			selectionTree=NULL;
		}
	}

	if(selectionTree==NULL){
		cout<<"selectionTree does not exists, close file"<<endl;
		delete selectionFile;
		cout<<"."<<endl;
		cout<<selectionfilepath.str().c_str()<<endl;
		selectionFile=new TFile(selectionfilepath.str().c_str(),"RECREATE");
		cout<<"."<<endl;
		this->selectionTree=new TTree("selectionTree",treeDescription.str().c_str());
		cout<<"."<<endl;
		createdNewTree=true;
		cout<<"\n\n\n***************************************************************\n";
		cout<<"***************************************************************\n";
		cout<<"there exists no tree:\'selectionTree\"\tcreate new one."<<selectionTree<<"\n";
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
	nDiamondHits=eventReader->getNClusters(8);
	isInFiducialCut=true;
	if(hasValidSiliconTrack){
		Float_t fiducialValueX=0;
		Float_t fiducialValueY=0;
		for(UInt_t plane=0;plane<4;plane++){
			fiducialValueX+=eventReader->getCluster(plane*2,0).getPosition();
			fiducialValueY+=eventReader->getCluster(plane*2+1,0).getPosition();
		}
		fiducialValueX/=4.;
		fiducialValueY/=4.;
		isInFiducialCut=isInFiducialCut&&fiducialValueX>settings->getSi_avg_fidcut_xlow();
		isInFiducialCut=isInFiducialCut&&fiducialValueX<settings->getSi_avg_fidcut_xhigh();
		isInFiducialCut=isInFiducialCut&&fiducialValueX>settings->getSi_avg_fidcut_ylow();
		isInFiducialCut=isInFiducialCut&&fiducialValueX<settings->getSi_avg_fidcut_yhigh();
//		cout<<"fidCut:"<<fiducialValueX<<"/"<<fiducialValueY<<":"<<isInFiducialCut<<endl;
	}
	else
		isInFiducialCut=false;
//		UInt_t nDiamondHits; //number of clusters in diamond plane;
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

void TSelectionClass::setBranchAdressess(){
	selectionTree->Branch("nDiamondHits",&nDiamondHits,"nDiamondHits/i");
	selectionTree->Branch("isInFiducialCut",&isInFiducialCut,"isInFiducialCut/O");
	selectionTree->Branch("isDetMasked",&isDetMasked,"isDetMasked/O");
	selectionTree->Branch("hasValidSiliconTrack",&hasValidSiliconTrack,"hasValidSiliconTrack/O");
	selectionTree->Branch("isDiaMasked",&this->isDiaMasked,"isDiaMasked");
}
