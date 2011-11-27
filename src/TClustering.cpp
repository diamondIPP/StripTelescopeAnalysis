/*
 * TClustering.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 *
 *      UPDATE TCLUSTER_REVISION IN HEADER FILE IF YOU MAKE CHANGES IN TCLUSTER!!!!!
 */

#include "../include/TClustering.hh"

TClustering::TClustering(int runNumber,int seedDetSigma,int hitDetSigma,int seedDiaSigma, int hitDiaSigma) {
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
	this->seedDetSigma=seedDetSigma;
	this->hitDetSigma=hitDetSigma;
	this->seedDiaSigma=seedDiaSigma;
	this->hitDiaSigma=hitDiaSigma;
	this->runNumber=runNumber;
	verbosity=0;
	this->maxDetAdcValue=255;
	this->maxDiaAdcValue=4095;
	pVecvecCluster=&vecvecCluster;

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
	vecvecCluster.resize(9);
	if(!createClusterTree(nEvents)) return;
	setBranchAdresses();
	cout<<"\n\n******************************************\n";
	cout<<    "**************Start Clustering...*********\n";
	cout<<"******************************************\n\n"<<endl;
	cout<< "\tSNRs for silicon: "<<seedDetSigma<<"/"<<hitDetSigma<<endl;
	cout<< "\tSNRs for diamond: "<<seedDiaSigma<<"/"<<hitDiaSigma<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){

		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		clusterEvent();
		clusterTree->Fill();
	}
}

void TClustering::clusterEvent()
{

	//Cluster Silicon Planes
	for(unsigned int det=0;det<9;det++){
		//clear vecCluster
		vecCluster[det].clear();

		//cluster Plane
		clusterPlane(det);

		//fill clusters in vecvecCluster
		if(vecvecCluster.size()>det){
			vecvecCluster.at(det).clear();
			for(unsigned int cl=0;cl<vecCluster[det].size();cl++)
				vecvecCluster.at(det).push_back(vecCluster[det].at(cl));
		}
		else
			cout<<"Something is going wrong vecvecCluster is to small....."<<endl;
		//hNumberOfSeeds[det]->Fill(numberOfSeeds);
	}


}

void TClustering::clusterPlane(int det){
	nClusters[det]=0;
	int maxChannels= (det==8)?N_DIA_CHANNELS:N_DET_CHANNELS;
	for(int ch=0;ch<maxChannels;ch++){
		Float_t sigma=eventReader->getPedestalSigma(det,ch);
		Float_t signal = (Float_t)eventReader->getAdcValue(det,ch)-eventReader->getPedestalMean(det,ch);
		if(verbosity>2&&nEvent==0&&det==0&&ch<20)cout<<nEvent<<" "<<det<<" "<<ch<<" "<<signal<<" "<<sigma<<" "<<flush;
		//if(det==8)cout<<nEvent<<" # "<<det<<" # "<<ch<<" "<<signal<<" "<<sigma<<" "<<endl;
		if(sigma==0){
			if(verbosity>1)cout<<nEvent<<" # "<<det<<" # "<<ch<<" sigma==0"<<endl;
			continue;
		}
		Float_t SNR=signal/sigma;
		if(SNR!=eventReader->getSignalInSigma(det,ch))cout<<"in the SNR there is something wrong...";
		if(verbosity>2&&nEvent==0&&det==0&&ch<20)cout<<SNR<<endl;

		if(det<8 && SNR>this->seedDetSigma){
			if(verbosity>3)cout<<"Found a Seed "<<nEvent<<eventReader->getCurrent_event() <<" "<<det<<" "<<ch<<" "<<signal<<" "<<SNR<<" "<<flush;
			ch=combineCluster(det,ch,this->maxDetAdcValue);
			if(verbosity>3)cout<<"new channel no.:"<<ch<<endl;
		}
		if(det==8 && SNR>this->seedDiaSigma){
			if(verbosity>2)cout<<"Found a DiaSeed "<<nEvent<<" "<<det<<" "<<ch<<" "<<signal<<" "<<SNR<<" "<<eventReader->getCurrent_event()<<flush;
			ch=combineCluster(det,ch,this->maxDiaAdcValue);
			if(verbosity>2)cout<<"new channel no.:"<<ch<<endl;
		}
	}
	if(verbosity>1&&det==8){
		cout<<"Clustered Plane "<<det<<" with "<<nClusters[det]<<" "<<vecCluster[det].size()<<"with size: ";
		for(unsigned int cl=0;cl<vecCluster[det].size();cl++)
			cout<<" "<<vecCluster[det].at(cl).size()<<"/"<<vecCluster[det].at(cl).getCharge();
		cout<<endl;
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
int TClustering::combineCluster(int det, int ch,int maxAdcValue){
	if((verbosity>2&&det==8)||verbosity>3)cout<<"combine Cluster...start:"<<ch<<" ";

	Float_t sigma=eventReader->getPedestalSigma(det,ch);
	Float_t signal = (Float_t)eventReader->getAdcValue(det,ch)-eventReader->getPedestalMean(det,ch);
	Float_t adcValueInSigma=signal/sigma;

	//create Cluster
	int seedSigma;
	int hitSigma;
	if (det==8){
		hitSigma=this->hitDiaSigma;
		seedSigma=this->seedDiaSigma;
	}
	else{
		hitSigma=this->hitDetSigma;
		seedSigma=this->seedDetSigma;
	}
	TCluster cluster(nEvent,seedSigma,hitSigma);

	//look for hit channels smaller than or equal  to the seed channel
	if(verbosity>2)cout<<cluster.size()<<" ";
	for(int currentCh=ch;adcValueInSigma>hitSigma&&currentCh>=0;currentCh--){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		UShort_t adcValue=eventReader->getAdcValue(det,currentCh);
		signal = (Float_t)adcValue-eventReader->getPedestalMean(det,currentCh);
		adcValueInSigma=signal/sigma;
		if(sigma!=0&&adcValueInSigma>hitSigma){
			cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue);//todo add saturated
		}
		else{
			if((verbosity>2&&det==8)||verbosity>3)cout<<" ["<<currentCh<<"/"<<signal<<"/"<<sigma<<"/"<<adcValueInSigma<<"] ";
			break;
		}
	}
	if((verbosity>2&&det==8)||verbosity>3)cout<<" ."<<cluster.size()<<". ";
	int currentCh;
	for(currentCh=ch+1;currentCh<N_DET_CHANNELS;currentCh++){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		UShort_t adcValue=eventReader->getAdcValue(det,currentCh);
		if(sigma==0)cout<<"$";
		signal = (Float_t)adcValue-eventReader->getPedestalMean(det,currentCh);
		adcValueInSigma=signal/sigma;
		if(sigma!=0&&adcValueInSigma>hitSigma&&sigma!=0){
			cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue);
		}
		else{
			if((verbosity>2&&det==8)||verbosity>3)cout<<" ["<<currentCh<<"/"<<signal<<"/"<<adcValueInSigma<<"] ";
			break;
		}
	}
	cluster.checkCluster();
	vecCluster[det].push_back(cluster);
	nClusters[det]++;
	if((verbosity>2&&det==8)||verbosity>3)cout<<"\tclusterSize: "<<cluster.size()<<endl;
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
			cout<<"tree has enough entries....check Rev"<<endl;
			clusterRev=-1;
			cout<<"#";
			clusterTree->SetBranchAddress("clusterRev",&clusterRev);
			cout<<"#";
			clusterTree->GetEvent(0);
			cout<<"#";
			cout<<"ClusterTree has revision: rev."<<clusterRev<<" current rev."<<TCLUSTER_REVISION<<endl;
			if(clusterRev==TCLUSTER_REVISION)
				return false;
			else{
				cout<<"ClusterTree has wrong revision: rev."<<clusterRev<<" instead of rev."<<TCLUSTER_REVISION<<endl;
				clusterTree->Delete();
				clusterTree=NULL;
			}
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
		cout<<"\n\n\n***************************************************************\n";
		cout<<"***************************************************************\n";
		cout<<"there exists no tree:\'clusterTree\"\tcreate new one."<<clusterTree<<"\n";
		cout<<"***************************************************************\n";
		cout<<"***************************************************************\n"<<endl;
	}

	return createdNewTree;
}


void TClustering::setBranchAdresses(){
	cout<<"set Branch adresses..."<<endl;

	clusterRev=TCLUSTER_REVISION;
	clusterTree->Branch("eventNumber",&nEvent,"eventNumber/i");
	clusterTree->Branch("runNumber",&runNumber,"runNumber/i");
	clusterTree->Branch("nClusters",&nClusters,"nClusters/i[9]");
	clusterTree->Branch("clusterRev",&clusterRev,"clusterRev/I");
	//clusterTree->Branch("vecvecChannel",&vecvecChannel[0])
	// example t1.Branch("tracks","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pTracks);
	clusterTree->Branch("clusters","std::vector<std::vector<TCluster> >",&pVecvecCluster);
}






