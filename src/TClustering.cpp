/*
 * TClustering.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 *
 *      UPDATE TCLUSTER_REVISION IN HEADER FILE IF YOU MAKE CHANGES IN TCLUSTER!!!!!
 */

#include "../include/TClustering.hh"

TClustering::TClustering(TSettings* set){//int runNumber,int seedDetSigma,int hitDetSigma,int seedDiaSigma, int hitDiaSigma) {
	cout<<"**********************************************************"<<endl;
	cout<<"*************TClustering::TClustering*********************"<<endl;
	cout<<"**********************************************************"<<endl;
	// TODO Auto-generated constructor stub

	settings=NULL;
	if(set==0){
		cerr<< "Settings ==0 , exit"<<endl;
		exit(-1);
	}
	setSettings(set);
	UInt_t runNumber = settings->getRunNumber();
	sys = gSystem;
	settings->goToPedestalTreeDir();
	eventReader=new TADCEventReader(settings->getPedestalTreeFilePath(),settings);//settings->getRunNumber());
	histSaver=new HistogrammSaver(settings);
	settings->goToClusterAnalysisDir();
	stringstream plotsPath;
	plotsPath<<sys->pwd()<<"/";
	histSaver->SetPlotsPath(plotsPath.str().c_str());
	histSaver->SetRunNumber(runNumber);
	settings->goToPedestalTreeDir();
	this->runNumber=runNumber;
	verbosity=settings->getVerbosity();
	createdTree=false;
	pEvent=0;//new TEvent();
	for(UInt_t det=0;det<9;det++){
		stringstream histName;
		histName<<"hEtaDistribution_"<<det;//<<TADCEventReader::getStringForPlane(det);
		hEtaDistribution[det]=new TH1F(histName.str().c_str(),histName.str().c_str(),1024,0,1);
	}
	nInvalidReadout=0;
}

TClustering::~TClustering() {
	clusterFile->cd();

    if(clusterTree!=NULL&&this->createdTree){
        clusterFile->GetListOfKeys()->Print();

		cout<<"Invalid readouts: "<<nInvalidReadout<<endl;
		if(verbosity)cout<<"CLOSING TREE"<<endl;
		if(verbosity)cout<<"pedestalTree"<<" "<<settings->getPedestalTreeFilePath()<<" "<<filepath.str().c_str()<<endl;
		clusterTree->AddFriend("pedestalTree",settings->getPedestalTreeFilePath().c_str());
		if(verbosity)cout<<"rawTree"<<" "<<settings->getRawTreeFilePath()<<" "<<rawFilePath.str().c_str()<<endl;
		clusterTree->AddFriend("rawTree",settings->getRawTreeFilePath().c_str());
		if(verbosity)cout<<"save clusterTree: "<<clusterTree->GetListOfFriends()->GetEntries()<<endl;
		clusterTree->Write();
		saveEtaCorrections();
	}
	delete clusterFile;
	delete eventReader;
	delete histSaver;
	settings->goToOutputDir();
}

void TClustering::setSettings(TSettings* settings){
	this->settings = settings;
}

void TClustering::ClusterEvents(UInt_t nEvents)
{
	if(settings==NULL) settings=new TSettings("");//todo anpassen
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
		cout<< "\tSNRs for silicon plane "<<det<<": "<<settings->getClusterSeedFactor(det,0)<<"/"<<settings->getClusterHitFactor(det,0)<<endl;
	cout<<endl;
	for(UInt_t det=TPlaneProperties::getDetDiamond();det< TPlaneProperties::getNDetectors();det++)
		cout<< "\tSNRs for diamond plane "<<det<<": "<<settings->getClusterSeedFactor(det,0)<<"/"<<settings->getClusterHitFactor(det,0)<<endl;
	UInt_t validEvents=0;
	for(nEvent=0;nEvent<nEvents;nEvent++){

		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->LoadEvent(nEvent);
		clusterEvent();
		addToEtaDistributions();
		clusterTree->Fill();

		if (pEvent->isValidSiliconEvent())
			validEvents++;
	}

	cout<<"nvalid Events: "<<validEvents<<" of "<<nEvents<<endl;
}

void TClustering::clusterEvent()
{

	//Cluster Planes
	for(unsigned int det=0;det<TPlaneProperties::getNDetectors();det++){
		//clear vecCluster
		vecCluster[det].clear();

		//cluster Plane
		clusterDetector(det);
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
	if(pEvent->hasInvalidReadout()){
		if(verbosity>4)cout<<nEvent<<": InvalidReadout"<<endl;
		nInvalidReadout++;
	}

}

void TClustering::clusterDetector(UInt_t det){
	if(det>=TPlaneProperties::getNDetectors()){
		cerr<<"det is bigger than NDetectors... EXIT"<<endl;
		exit(-1);
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
			if(verbosity>8 ||(det ==8 && verbosity>3))cout<<nEvent<<" # "<<det<<" # "<<ch<<" sigma==0"<<endl;
			continue;
		}
		Float_t SNR=eventReader->getSignalInSigma(det,ch);
		if(SNR!=eventReader->getSignalInSigma(det,ch))cout<<"in the SNR there is something wrong...";
		//if(verbosity>2&&nEvent==0&&det==8&&ch<TPlaneProperties::getNChannels(det))cout<<SNR<<flush;


		if( SNR>settings->getClusterSeedFactor(det,ch)){
			if(verbosity>8||(det ==8 && verbosity>3))cout<<"Found a Seed "<<nEvent<<" "<<eventReader->getCurrent_event() <<" "<<det<<" "<<ch<<" "<<signal<<" "<<SNR<<" "<<flush;
			ch=combineCluster(det,ch);
			if(verbosity>20||(det ==8 && verbosity>5))cout<<"new channel no.:"<<ch<<flush;
		}
	}
	if(verbosity>8||(det ==8 && verbosity>3))cout<<endl;

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
int TClustering::combineCluster(UInt_t det, UInt_t ch){
//	int maxAdcValue = TPlaneProperties::getMaxSignalHeight(det); // DA
	int maxAdcValue = settings->getMaxSignalHeight(det); // DA
	if((verbosity>10&&det==8)||verbosity>11)cout<<"combine Cluster...start:"<<ch<<" ";

	Float_t sigma=eventReader->getPedestalSigma(det,ch);
	Float_t signal =eventReader->getSignal(det,ch);
	Float_t adcValueInSigma=eventReader->getSignalInSigma(det,ch);
	Int_t adcValue= eventReader->getAdcValue(det,ch);
	Float_t cmNoise = eventReader->getCMNoise(det,ch);

	//create Cluster
	int seedSigma=settings->getClusterSeedFactor(det,ch);
	int hitSigma=settings->getClusterHitFactor(det,ch);
	bool isScreened=false;
	int maxChannel=TPlaneProperties::getNChannels(det);
	TCluster cluster(nEvent,(UChar_t)det,seedSigma,hitSigma,maxChannel,cmNoise);

	//look for hit channels smaller than or equal  to the seed channel
	if(verbosity>10)cout<<cluster.size()<<" ";
	UInt_t currentCh;
	for(currentCh=ch;adcValueInSigma>hitSigma&&currentCh>=0&&currentCh<=TPlaneProperties::getNChannels(det);currentCh--){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		adcValue=eventReader->getAdcValue(det,currentCh);
		if(verbosity&&sigma<=0)cout<<currentCh<<":sigma<0 ";
		signal =eventReader->getSignal(det,currentCh);
		adcValueInSigma=eventReader->getSignalInSigma(det,currentCh);
		isScreened=this->settings->isDet_channel_screened(det,currentCh)||adcValue==maxAdcValue;
		if(sigma!=0&&adcValueInSigma>hitSigma){
			float pedMean = eventReader->getPedestalMean(det,currentCh,false);
			float pedMeanCMN = eventReader->getPedestalMean(det,currentCh,true);
			float pedSigma = eventReader->getPedestalSigma(det,currentCh,false);
			float pedSigmaCMN = eventReader->getPedestalSigma(det,currentCh,true);
//			cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,TPlaneProperties::isSaturated(det,adcValue),isScreened); // DA
			cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,settings->isSaturated(det,adcValue),isScreened); // DA
			//			cluster.addChannel(currentCh,signal,adcValueInSigma,adcValue,adcValue>=maxAdcValue,isScreened);//todo add saturated
		}
		else{
			if((verbosity>10&&det==8)||verbosity>11)cout<<" ["<<currentCh<<"/"<<signal<<"/"<<sigma<<"/"<<adcValueInSigma<<"] ";
			break;
		}
	}
	if(currentCh>=0){
		float pedMean = eventReader->getPedestalMean(det,currentCh,false);
		float pedMeanCMN = eventReader->getPedestalMean(det,currentCh,true);
		float pedSigma = eventReader->getPedestalSigma(det,currentCh,false);
		float pedSigmaCMN = eventReader->getPedestalSigma(det,currentCh,true);
		if(currentCh>=0&&currentCh<TPlaneProperties::getNChannels(det))
//			cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,TPlaneProperties::isSaturated(det,adcValue),isScreened); // DA
			cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,settings->isSaturated(det,adcValue),isScreened); // DA
	}
	if((verbosity>10&&det==8)||verbosity>11)cout<<" ."<<cluster.size()<<". ";
	for(currentCh=ch+1;currentCh<TPlaneProperties::getNChannels(det);currentCh++){
		sigma=eventReader->getPedestalSigma(det,currentCh);
		adcValue=eventReader->getAdcValue(det,currentCh);
		if(verbosity&&sigma<=0)cout<<currentCh<<":sigma<0 ";
		signal =eventReader->getSignal(det,currentCh);
		adcValueInSigma=eventReader->getSignalInSigma(det,currentCh);
		isScreened=this->settings->isDet_channel_screened(det,currentCh);
		if(sigma!=0&&adcValueInSigma>hitSigma&&sigma!=0){
			float pedMean = eventReader->getPedestalMean(det,currentCh,false);
			float pedMeanCMN = eventReader->getPedestalMean(det,currentCh,true);
			float pedSigma = eventReader->getPedestalSigma(det,currentCh,false);
			float pedSigmaCMN = eventReader->getPedestalSigma(det,currentCh,true);
//			cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,TPlaneProperties::isSaturated(det,adcValue),isScreened); // DA
			cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,settings->isSaturated(det,adcValue),isScreened); // DA
		}
		else{
			if((verbosity>10&&det==8)||verbosity>11)cout<<" ["<<currentCh<<"/"<<signal<<"/"<<adcValueInSigma<<"] ";
			break;
		}
	}
	if(currentCh<TPlaneProperties::getNChannels(det)){
		float pedMean = eventReader->getPedestalMean(det,currentCh,false);
		float pedMeanCMN = eventReader->getPedestalMean(det,currentCh,true);
		float pedSigma = eventReader->getPedestalSigma(det,currentCh,false);
		float pedSigmaCMN = eventReader->getPedestalSigma(det,currentCh,true);
//		cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,TPlaneProperties::isSaturated(det,adcValue),isScreened); // DA
		cluster.addChannel(currentCh,pedMean,pedSigma,pedMeanCMN,pedSigmaCMN,adcValue,settings->isSaturated(det,adcValue),isScreened); // DA
	}
	cluster.checkCluster();
	vecCluster[det].push_back(cluster);
	nClusters[det]++;
	if((verbosity>10&&det==8)||verbosity>11)cout<<"\tclusterSize: "<<cluster.size()<<endl;
	if(verbosity>8||(det ==8 && verbosity>3))cluster.Print();
	return currentCh;
}

bool TClustering::createClusterTree(int nEvents)
{
	bool createdNewFile=false;
	bool createdNewTree=false;
	stringstream clusterfilepath;
	clusterfilepath<<sys->pwd();
	clusterfilepath<<"/clusterData."<<runNumber<<".root";
	if(verbosity)cout<<"Try to open \""<<clusterfilepath.str()<<"\""<<endl;
	clusterFile=new TFile(clusterfilepath.str().c_str(),"READ");
	if(clusterFile->IsZombie()){
		delete clusterFile;
		if(verbosity)cout<<"clusterfile does not exist, create new one..."<<endl;
		createdNewFile =true;
		clusterFile= new TFile(clusterfilepath.str().c_str(),"RECREATE");
		clusterFile->cd();
	}
	else{
		createdNewFile=false;
		if(verbosity)cout<<"File exists"<<endl;
	}
	clusterFile->cd();
	stringstream treeDescription;
	treeDescription<<"Cluster Data of run "<<runNumber;
	if(!createdNewFile)
		clusterFile->GetObject("clusterTree",clusterTree);
	else
		clusterTree=NULL;
	if(clusterTree!=NULL){
		if(verbosity)cout<<"File and Tree Exists... \t"<<flush;
		if(clusterTree->GetEntries()>=nEvents){
			createdNewTree=false;
			if(verbosity)cout<<"tree has enough entries....check Rev"<<endl;
			clusterRev=-1;
			if(verbosity)cout<<"#";
			clusterTree->SetBranchAddress("clusterRev",&clusterRev);
			if(verbosity)cout<<"#";
			clusterTree->GetEvent(0);
			if(verbosity)cout<<"#";
			if(verbosity)cout<<"ClusterTree has revision: rev."<<clusterRev<<" current rev."<<TCluster::TCLUSTER_REVISION()<<endl;
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
		if(!createdNewFile||!clusterFile){
			clusterFile->Close();
			clusterFile=new TFile(clusterfilepath.str().c_str(),"RECREATE");
		}
		clusterFile->cd();
		this->clusterTree=new TTree("clusterTree",treeDescription.str().c_str());
		createdNewTree=true;
		if(verbosity)cout<<"\n\n***************************************************************\n";
		if(verbosity)cout<<"there exists no tree:\'clusterTree\"\tcreate new one."<<clusterTree<<"\n";
		if(verbosity)cout<<"***************************************************************\n"<<endl;
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
		TH1F *histo= createEtaIntegral(hEtaDistribution[det],histName.str());
		file->cd();
		histo->Write();
		hEtaDistribution[det]->Write();
	}
	file->Close();
}

void TClustering::setBranchAdresses(){
	if(verbosity)cout<<"set Branch adresses..."<<endl;

	clusterRev=TCluster::TCLUSTER_REVISION();
	if(verbosity)cout<<"Branch eventNumber"<<endl;
	clusterTree->Branch("eventNumber",&nEvent,"eventNumber/i");
	if(verbosity)cout<<"Branch runNumber"<<endl;
	clusterTree->Branch("runNumber",&runNumber,"runNumber/i");
	if(verbosity)cout<<"Branch nClusters"<<endl;
	clusterTree->Branch("nClusters",&nClusters,"nClusters/i[9]");
	if(verbosity)cout<<"Branch clusterRev"<<endl;
	clusterTree->Branch("clusterRev",&clusterRev,"clusterRev/i");
	if(verbosity)cout<<"Branch clusters"<<endl;
	//clusterTree->Branch("vecvecChannel",&vecvecChannel[0])
	// example t1.Branch("tracks","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pTracks);
	//	clusterTree->Branch("clusters","std::vector<std::vector<TCluster> >",&pVecvecCluster);
	if(verbosity)cout<<"Branch event"<<endl;
	pEvent=0;
	clusterTree->Branch("event","TEvent",&pEvent);
}


TH1F *TClustering::createEtaIntegral(TH1F *histo, std::string histName)
{
	UInt_t nBins = histo->GetNbinsX();
	TH1F *hIntegral=new TH1F(histName.c_str(),histName.c_str(),nBins,0,1);
	Int_t entries = histo->GetEntries();
	entries -=  histo->GetBinContent(0);
	entries -=  histo->GetBinContent(nBins+1);
	Int_t sum =0;
	for(UInt_t bin=1;bin<nBins+1;bin++){
		Int_t binContent = histo->GetBinContent(bin);
		sum +=binContent;
		Float_t pos =  histo->GetBinCenter(bin);
		hIntegral->Fill(pos, (Float_t)sum/(Float_t)entries);
	}
	return hIntegral;
}


