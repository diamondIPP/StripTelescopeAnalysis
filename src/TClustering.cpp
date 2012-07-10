/*
 * TClustering.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 *
 *      UPDATE TCLUSTER_REVISION IN HEADER FILE IF YOU MAKE CHANGES IN TCLUSTER!!!!!
 */

#include "../include/TClustering.hh"

TClustering::TClustering(TSettings* settings){
	cout<<"**********************************************************"<<endl;
	cout<<"*************TClustering::TClustering*********************"<<endl;
	cout<<"**********************************************************"<<endl;
	if(settings==0)
		settings=new TSettings();
	setSettings(settings);
	UInt_t runNumber = settings->getRunNumber();
	sys = gSystem;

  sys->MakeDirectory(settings->getRelativePath().c_str());
  sys->cd(settings->getRelativePath().c_str());

	rawFilePath<<"rawData."<<runNumber<<".root";
	filepath.str("");
	filepath<<"pedestalData."<<runNumber<<".root";
	cout<<"currentPath: "<<sys->pwd()<<endl;
	cout<<filepath.str()<<endl;
	eventReader=new TADCEventReader(filepath.str(),settings->getRunNumber());
	histSaver=new HistogrammSaver();
	sys->MakeDirectory("clustering");
	sys->cd("clustering");
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	sys->cd("..");
	this->runNumber=runNumber;
	verbosity=0;
	this->maxDetAdcValue=255;
	this->maxDiaAdcValue=4095;
//	pVecvecCluster=&vecvecCluster;
	settings=NULL;
	createdTree=false;
	pEvent=0;//new TEvent();
	for(UInt_t det=0;det<9;det++){
		stringstream histName;
		histName<<"hEtaDistribution_"<<det;//<<TADCEventReader::getStringForPlane(det);
		hEtaDistribution[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,0,1);
	}
}

TClustering::~TClustering() {
	clusterFile->cd();
	if(clusterTree!=NULL&&this->createdTree){
		cout<<"CLOSING TREE"<<endl;
		cout<<eventReader->getTree()->GetName()<<" "<<filepath.str().c_str()<<endl;
		clusterTree->AddFriend(eventReader->getTree()->GetName(),filepath.str().c_str());
		cout<<"rawTree"<<" "<<rawFilePath.str().c_str()<<endl;
		clusterTree->AddFriend("rawTree",rawFilePath.str().c_str());
		cout<<"save clusterTree: "<<clusterTree->GetListOfFriends()->GetEntries()<<endl;
		clusterTree->Write();
		saveEtaCorrections();
	}
	//clusterTree->Delete();
	delete clusterFile;
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TClustering::setSettings(TSettings* settings){
	this->settings = settings;
}

void TClustering::ClusterEvents(UInt_t nEvents){
	if(settings==NULL) settings=new TSettings("");
//	vecvecCluster.resize(9);
	createdTree=createClusterTree(nEvents);
	if(!createdTree) return;
	setBranchAdresses();
	cout<<"\n\n******************************************\n";
	cout<<    "**************Start Clustering...*********\n";
	cout<<"******************************************\n\n"<<endl;
	for(UInt_t det=0;det<9;det++){
		if(hEtaDistribution[det]!=0)delete hEtaDistribution[det];
		stringstream histName;
		histName<<"hEtaDistribution_"<<det;//<<TADCEventReader::getStringForPlane(det);
		hEtaDistribution[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),nEvents/50,0,1);
	}
	for(UInt_t det=0;det< TPlaneProperties::getNSiliconDetectors();det++)
		cout<< "\tSNRs for silicon plane "<<det<<": "<<settings->getClusterSeedFactor(det)<<"/"<<settings->getClusterHitFactor(det)<<endl;
	cout<<endl;
	for(UInt_t det=TPlaneProperties::getDetDiamond();det< TPlaneProperties::getNDetectors();det++)
			cout<< "\tSNRs for diamond plane "<<det<<": "<<settings->getClusterSeedFactor(det)<<"/"<<settings->getClusterHitFactor(det)<<endl;
	UInt_t validEvents=0;
	for(nEvent=0;nEvent<nEvents;nEvent++){

		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
		clusterEvent();
		addToEtaDistributions();
		clusterTree->Fill();

		if(pEvent->isValidSiliconEvent())
			validEvents++;
	}

	cout<<"'nvalid Events: "<<validEvents<<" of "<<nEvents<<endl;
}

void TClustering::clusterEvent()
{

	//Cluster Planes
	for(unsigned int det=0;det<TPlaneProperties::getNDetectors();det++){
		//clear vecCluster
		vecCluster[det].clear();

		//cluster Plane
		clusterDetector(det);
		//fill clusters in vecvecCluster
//		if(verbosity>10)cout<<"fill "<<det<<" vecvecCluster "<<flush;
//		if(det<vecvecCluster.size()){
//			if(verbosity>10)cout<<vecvecCluster.size()<<"."<<flush;
//			vecvecCluster.at(det).clear();
//			if(verbosity>10)cout<<","<<vecCluster[det].size()<<flush;
//			for(unsigned int cl=0;cl<vecCluster[det].size();cl++){
//				if(verbosity>10)cout<<"."<<cl<<flush;
//				TCluster cluster=vecCluster[det].at(cl);
//				if(verbosity>10)cout<<"."<<cl<<flush;
//				vecvecCluster.at(det).push_back(cluster);
//			}
//		}
//		else
//			cout<<"Something is going wrong vecvecCluster is to small....."<<endl;
//		//hNumberOfSeeds[det]->Fill(numberOfSeeds);
//		if(verbosity>10)cout<<"Done with detector "<<det<<endl;
	}

	//Save Planes to Event
	if(pEvent!=NULL) {delete pEvent;pEvent=NULL;}
	pEvent = new TEvent(nEvent);
	if(verbosity>10)cout<<"."<<flush;

	//siliconPlanes
	for(UInt_t nplane=0;nplane<TPlaneProperties::getNSiliconPlanes();nplane++){
		TPlane plane(nplane,vecCluster[nplane*2],vecCluster[nplane*2+1],TPlaneProperties::kSilicon);
		if(verbosity>10)plane.Print(1);
		pEvent->addPlane(plane,nplane);
		if(verbosity>10)cout<<nplane<<"."<<flush;
	}

	//diamondPlanes
	TPlane plane(TPlaneProperties::getDiamondPlane(),vecCluster[TPlaneProperties::getDetDiamond()],TPlaneProperties::kDiamond);
	if(verbosity>10)cout<<4<<"."<<flush;
	pEvent->addPlane(plane,TPlaneProperties::getDiamondPlane());
	if(true){pEvent->isValidSiliconEvent();}
	if(verbosity>8){
		cout<<"\n"<<nEvent<<" "<<pEvent->getEventNumber()<<" "<<pEvent->isValidSiliconEvent()<<" ";
//		for (UInt_t det=0;det<vecvecCluster.size();det++){
//			cout<<vecvecCluster.at(det).size()<<" ";
//		}
		cout<<endl;
	}

}

void TClustering::clusterDetector(int det){
	if(det>=TPlaneProperties::getNDetectors()){
		cerr<<"det is bigger than NDetectors... EXIT"<<endl;
		exit;
	}
	nClusters[det]=0;
	int maxChannels= TPlaneProperties::getNChannels(det);
	if(verbosity>10)cout<<"ClusterDetector"<<det<<" "<<maxChannels<<endl;
	for(int ch=0;ch<maxChannels;ch++){
		//if(verbosity>30&&nEvent==0&&det==8&&ch<128)cout<<nEvent<<flush;

		Float_t sigma=eventReader->getPedestalSigma(det,ch);
		Float_t signal = eventReader->getSignal(det,ch);

		//if(verbosity>9&&nEvent==0&&det==8&&ch<128)cout<<" "<<det<<" "<<ch<<" "<<signal<<" "<<sigma<<" "<<flush;
		//if(det==8)cout<<nEvent<<" # "<<det<<" # "<<ch<<" "<<signal<<" "<<sigma<<" "<<endl;
		if(sigma==0){
			if(verbosity>1)cout<<nEvent<<" # "<<det<<" # "<<ch<<" sigma==0"<<endl;
			continue;
		}
		Float_t SNR=eventReader->getSignalInSigma(det,ch);
		if(SNR!=eventReader->getSignalInSigma(det,ch))cout<<"in the SNR there is something wrong...";
		//if(verbosity>2&&nEvent==0&&det==8&&ch<TPlaneProperties::getNChannels(det))cout<<SNR<<flush;


		if( SNR>settings->getClusterSeedFactor(det)){
			if(verbosity>3)cout<<"Found a Seed "<<nEvent<<" "<<eventReader->getCurrent_event() <<" "<<det<<" "<<ch<<" "<<signal<<" "<<SNR<<" "<<flush;
			ch=combineCluster(det,ch,TPlaneProperties::getMaxSignalHeight(det));
			if(verbosity>20)cout<<"new channel no.:"<<ch<<flush;
		}
	}
	if(verbosity>3)cout<<endl;

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
	if((verbosity>10&&det==8)||verbosity>11)cout<<"combine Cluster...start:"<<ch<<" ";

	Float_t sigma=eventReader->getPedestalSigma(det,ch);
	Float_t signal =eventReader->getSignal(det,ch);
	Float_t adcValueInSigma=eventReader->getSignalInSigma(det,ch);
	UShort_t adcValue=(UShort_t)eventReader->getAdcValue(det,ch);

	//create Cluster
	int seedSigma=settings->getClusterSeedFactor(det);
	int hitSigma=settings->getClusterHitFactor(det);
	bool isScreened;
	int maxChannel=TPlaneProperties::getNChannels(det);

	TCluster cluster(nEvent,(UChar_t)det,seedSigma,hitSigma,maxChannel);

	//look for hit channels smaller than or equal  to the seed channel
	if(verbosity>10)cout<<cluster.size()<<" ";
	UInt_t currentCh;
	for(currentCh=ch;adcValueInSigma>hitSigma&&currentCh>=0;currentCh--){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		adcValue=eventReader->getAdcValue(det,currentCh);
		if(verbosity&&sigma<=0)cout<<currentCh<<":sigma<0 ";
		signal =eventReader->getSignal(det,currentCh);
		adcValueInSigma=eventReader->getSignalInSigma(det,currentCh);
		isScreened=this->settings->isDet_channel_screened(det,currentCh)||adcValue==maxAdcValue;
		if(sigma!=0&&adcValueInSigma>hitSigma){
			cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue,isScreened);//todo add saturated
		}
		else{
			if((verbosity>10&&det==8)||verbosity>11)cout<<" ["<<currentCh<<"/"<<signal<<"/"<<sigma<<"/"<<adcValueInSigma<<"] ";
			break;
		}
	}
	if(currentCh>=0)
		cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue,isScreened);//todo add saturated
	if((verbosity>10&&det==8)||verbosity>11)cout<<" ."<<cluster.size()<<". ";
	for(currentCh=ch+1;currentCh<TPlaneProperties::getNChannels(det);currentCh++){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		adcValue=eventReader->getAdcValue(det,currentCh);
		if(verbosity&&sigma<=0)cout<<currentCh<<":sigma<0 ";
		signal =eventReader->getSignal(det,currentCh);
		adcValueInSigma=eventReader->getSignalInSigma(det,currentCh);
		isScreened=this->settings->isDet_channel_screened(det,currentCh);
		if(sigma!=0&&adcValueInSigma>hitSigma&&sigma!=0){
			cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue,isScreened);
		}
		else{
			if((verbosity>10&&det==8)||verbosity>11)cout<<" ["<<currentCh<<"/"<<signal<<"/"<<adcValueInSigma<<"] ";
			break;
		}
	}
	if(currentCh<TPlaneProperties::getNChannels(det)){
		cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue,isScreened);//todo add saturated
	}
	cluster.checkCluster();
	vecCluster[det].push_back(cluster);
	nClusters[det]++;
	if((verbosity>10&&det==8)||verbosity>11)cout<<"\tclusterSize: "<<cluster.size()<<endl;
	if(verbosity>1)cluster.Print();
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
		delete clusterFile;
		cout<<"clusterfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		clusterFile= new TFile(clusterfilepath.str().c_str(),"RECREATE");
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
			cout<<"ClusterTree has revision: rev."<<clusterRev<<" current rev."<<TCluster::TCLUSTER_REVISION()<<endl;
			if(clusterRev==TCluster::TCLUSTER_REVISION())
				return false;
			else{
				cout<<"ClusterTree has wrong revision: rev."<<clusterRev<<" instead of rev."<<TCluster::TCLUSTER_REVISION()<<endl;
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
		cout<<"\n\n***************************************************************\n";
		cout<<"there exists no tree:\'clusterTree\"\tcreate new one."<<clusterTree<<"\n";
		cout<<"***************************************************************\n"<<endl;
	}

	return createdNewTree;
}

void TClustering::addToEtaDistributions()
{
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		for(UInt_t cl=0;cl<pEvent->getNClusters(det);cl++){
			Float_t eta = pEvent->getCluster(det,cl).getEta();
			if(eta>=0&&eta<=1)
				hEtaDistribution[det]->Fill(eta);
		}
}

void TClustering::saveEtaCorrections(){
	stringstream etaCorFileName;
	etaCorFileName<<"etaCorrection."<<settings->getRunNumber()<<".root";
	TFile* file = new TFile(etaCorFileName.str().c_str(),"RECREATE");
	file->cd();
	for(UInt_t det=0;det<9;det++){
		stringstream histName;
		histName<<"hEtaIntegral_"<<det;
		UInt_t nBins = hEtaDistribution[det]->GetNbinsX();
		TH1F *histo=new TH1F(histName.str().c_str(),histName.str().c_str(),nBins,0,1);
		Int_t entries = hEtaDistribution[det]->GetEntries();
		entries -=  hEtaDistribution[det]->GetBinContent(0);
		entries -=  hEtaDistribution[det]->GetBinContent(nBins+1);
		Int_t sum =0;
		for(UInt_t bin=1;bin<nBins+1;bin++){
			Int_t binContent = hEtaDistribution[det]->GetBinContent(bin);
			sum +=binContent;
			Float_t pos =  hEtaDistribution[det]->GetBinCenter(bin);
			histo->Fill(pos, (Float_t)sum/(Float_t)entries);
		}
		file->cd();
		histo->Write();
		hEtaDistribution[det]->Write();
	}
	file->Close();
}

void TClustering::setBranchAdresses(){
	cout<<"set Branch adresses..."<<endl;

	clusterRev=TCluster::TCLUSTER_REVISION();
	cout<<"Branch eventNumber"<<endl;
	clusterTree->Branch("eventNumber",&nEvent,"eventNumber/i");
	cout<<"Branch runNumber"<<endl;
	clusterTree->Branch("runNumber",&runNumber,"runNumber/i");
	cout<<"Branch nClusters"<<endl;
	clusterTree->Branch("nClusters",&nClusters,"nClusters/i[9]");
	cout<<"Branch clusterRev"<<endl;
	clusterTree->Branch("clusterRev",&clusterRev,"clusterRev/i");
	cout<<"Branch clusters"<<endl;
	//clusterTree->Branch("vecvecChannel",&vecvecChannel[0])
	// example t1.Branch("tracks","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pTracks);
//	clusterTree->Branch("clusters","std::vector<std::vector<TCluster> >",&pVecvecCluster);
	cout<<"Branch event"<<endl;
	pEvent=0;
	clusterTree->Branch("event","TEvent",&pEvent);
}






