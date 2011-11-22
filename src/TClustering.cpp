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
	stringstream  filepath;
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
	verbosity=0;
}

TClustering::~TClustering() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TClustering::ClusterEvents(int nEvents)
{
	cout<<"\n\n******************************************\n";
	cout<<    "**************Start Clustering...*********\n";
	cout<<"******************************************\n\n"<<endl;
	for(nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		ClusterEvent();
	}
}

void TClustering::ClusterEvent()
{
	vecCluster.clear();
	for(int det=0;det<8;det++){
		//clear vecCluster
		for(int i=0;i<vecCluster.size();i++)
			delete vecCluster.at(i);
		vecCluster.clear();

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
		if(verbosity)cout<<"found "<<vecCluster.size()<<" Clusters in plane "<<det<<endl;
		//hNumberOfSeeds[det]->Fill(numberOfSeeds);
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
	if(verbosity>2)cout<<"\tclusterSize: "<<cluster->size()<<endl;
	return currentCh;
}





