/*
 * TClustering.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#include "../include/TClustering.hh"

TClustering::TClustering(int runNumber,int seedSigma,int hitSigma) {
	// TODO Auto-generated constructor stub
	sys = gSystem;
	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());

	sys->cd(runString.str().c_str());
	filepath.str("");
	filepath<<"pedestalData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str());
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("deadChannels");
	sys->cd("deadChannels");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
	this->seedSigma=seedSigma;
	this->hitSigma=hitSigma;
	this->runNumber=runNumber;
	verbosity=4;
}

TClustering::~TClustering() {
	// TODO Auto-generated destructor stub
	clusterFile->cd();
	clusterTree->Write();
	if(clusterTree!=NULL){
		clusterTree->AddFriend(eventReader->getTree()->GetName(),filepath.str().c_str());
		cout<<"save clusterTree: "<<clusterTree->GetListOfFriends()->GetEntries()<<endl;
	}
	clusterTree->Delete();
	clusterFile->Close();
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TClustering::ClusterEvents(int nEvents)
{
	if(!createClusterTree(nEvents)) return;
	setBranchAdresses();
	vecvecCluster.resize(9);
	cout<<"\n\n******************************************\n";
	cout<<    "**************Start Clustering...*********\n";
	cout<<"******************************************\n\n"<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){

		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		clusterEvent();
		clusterTree->Fill();
	}
}

void TClustering::clusterEvent()
{
	vecCluster.clear();
	for(int det=0;det<8;det++){
		if(verbosity)cout<<"."<<flush;
		//clear vecCluster
		for(int i=0;i<vecCluster.size();i++){
			if(verbosity)cout<<"."<<flush;
			delete vecCluster.at(i);
		}
		if(verbosity)cout<<"#"<<flush;
		for(int i=0;i<vecvecCluster.at(det).size();i++){
			delete vecvecCluster.at(det).at(i);
			if(verbosity)cout<<"/"<<flush;
		}
		if(verbosity)cout<<"."<<flush;
		vecvecCluster.at(det).clear();
		vecCluster.clear();
		if(verbosity)cout<<"."<<flush;

		clusterPlane(det);
		nClusters[det]=(UInt_t)vecCluster.size();
		if(verbosity)cout<<nEvent<<" found "<<vecCluster.size()<<" "<<vecvecCluster.at(det).size()<<" Clusters in plane "<<det<<endl;
		//hNumberOfSeeds[det]->Fill(numberOfSeeds);
	}

}

void TClustering::clusterPlane(int det){

	for(int ch=0;ch<N_DET_CHANNELS;ch++){
		Float_t sigma=eventReader->getPedestalSigma(det,ch);
		Float_t signal = (Float_t)eventReader->getDet_ADC(det,ch)-eventReader->getPedestalMean(det,ch);
		if(verbosity>2&&nEvent==0&&det==0&&ch<20)cout<<nEvent<<" "<<det<<" "<<ch<<" "<<signal<<" "<<sigma<<" "<<flush;

		if(sigma==0){
			if(verbosity>1)cout<<nEvent<<" "<<det<<" "<<ch<<" sigma==0"<<endl;
			continue;
		}
		Float_t adcValueInSigma=signal/sigma;
		if(verbosity>2&&nEvent==0&&det==0&&ch<20)cout<<adcValueInSigma<<endl;

		if(adcValueInSigma>this->seedSigma){
			if(verbosity>2)cout<<"Found a Seed "<<nEvent<<" "<<det<<" "<<ch<<" "<<signal<<" "<<adcValueInSigma<<" "<<eventReader->getCurrent_event()<<flush;
			ch=combineCluster(det,ch);
			if(verbosity>2)cout<<"new channel no.:"<<ch<<endl;
		}
	}
}

/**
 * \brief combines all channels aorund channel ch which are higher than hitSigma to one Cluster,
 *
 * 			This function gets a position of a seed in terms of detector and channel
 * 			It combines all channels which are higher than hitSigma to one Cluster
 * 			first a cluster is created, than
 * 			in the first loop it the cluster is filled with all channels which are smaller
 * 			than the seed channel but have an signal to noise ratio higher than hitSigma.
 * 			If one channel has a SNR smaller than hitSigma the loop stops
 * 			In the second loop the function looks for channels which are higher than the
 * 			seed channel and have a signal in sigma bigger than
 *
 * 			\param det current detector where to combin the cluster
 * 			\param ch  channel which is seed and triggered the cluster
 *
 * 			\return first channel which is not part of the cluster
 */
int TClustering::combineCluster(int det, int ch){
	if(verbosity>2)cout<<"combine Cluster...";
	Float_t sigma=eventReader->getPedestalSigma(det,ch);
	Float_t signal = (Float_t)eventReader->getDet_ADC(det,ch)-eventReader->getPedestalMean(det,ch);
	Float_t adcValueInSigma=signal/sigma;

	//create Cluster
	TCluster *cluster =new TCluster(nEvent,this->seedSigma,this->hitSigma);

	//look for hit channels smaller than or equal  to the seed channel
	if(verbosity>2)cout<<cluster->size()<<" ";
	for(int currentCh=ch;adcValueInSigma>hitSigma&&currentCh>=0;currentCh--){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		signal = (Float_t)eventReader->getDet_ADC(det,currentCh)-eventReader->getPedestalMean(det,currentCh);
		adcValueInSigma=signal/sigma;
		if(sigma!=0&&adcValueInSigma>hitSigma)cluster->addChannel(currentCh,signal,adcValueInSigma);
	}
	if(verbosity>2)cout<<cluster->size()<<" ";
	int currentCh;
	for(currentCh=ch+1;currentCh<N_DET_CHANNELS&&adcValueInSigma>hitSigma;currentCh++){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		signal = (Float_t)eventReader->getDet_ADC(det,currentCh)-eventReader->getPedestalMean(det,currentCh);
		adcValueInSigma=signal/sigma;
		if(sigma!=0&&adcValueInSigma>hitSigma)cluster->addChannel(currentCh,signal,adcValueInSigma);
	}
	if(verbosity>2)cout<<cluster->size()<<" "<<currentCh<<" "<<ch<<" ";
	vecCluster.push_back(cluster);
	vecvecCluster.at(det).push_back((TCluster*)cluster->Clone());
	if(verbosity>2)cout<<"\tclusterSize: "<<cluster->size()<<endl;
	return currentCh;
}

bool TClustering::createClusterTree(int nEvents)
{
	bool createdNewFile=false;
	bool createdNewTree=false;
	stringstream clusterfilepath;
	clusterfilepath<<sys->pwd();
	clusterfilepath<<"/clusterData."<<runNumber<<".root";
	cout<<"Try to open \""<<clusterfilepath.str()<<"\""<<endl;
	clusterFile=new TFile(clusterfilepath.str().c_str(),"READ");
	if(clusterFile->IsZombie()){
		cout<<"clusterfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		clusterFile= new TFile(clusterfilepath.str().c_str(),"CREATE");
		clusterFile->cd();
	}
	else{
		createdNewFile=false;
		cout<<"File exists"<<endl;
	}
	clusterFile->cd();
	stringstream treeDescription;
	treeDescription<<"Cluster Data of run "<<runNumber;
	clusterFile->GetObject("clusterTree",clusterTree);
	if(clusterTree!=NULL){
		cout<<"File and Tree Exists... \t"<<flush;
		if(clusterTree->GetEntries()>=nEvents){
			createdNewTree=false;
			cout<<"tree has enough entries...."<<endl;
			return false;
		}
		else{
			clusterTree->Delete();
			clusterTree=NULL;
		}
	}
	if(clusterTree==NULL){
		clusterFile->Close();
		clusterFile=new TFile(clusterfilepath.str().c_str(),"RECREATE");
		this->clusterTree=new TTree("clusterTree",treeDescription.str().c_str());
		createdNewTree=true;
		cout<<"there exists no tree:\'clusterTree\"\tcreate new one."<<clusterTree<<endl;
	}

	return createdNewTree;
}


void TClustering::setBranchAdresses(){
	cout<<"set Branch adresses..."<<endl;
	clusterTree->Branch("eventNumber",&nEvent,"eventNumber/i");
	clusterTree->Branch("runNumber",&runNumber,"runNumber/i");
	clusterTree->Branch("nClusters",&nClusters,"nClusters/i[9]");
	TBranch *clusterBranch = clusterTree->Branch("clusters",&this->vecvecCluster);
}






