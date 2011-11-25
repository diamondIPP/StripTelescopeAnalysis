/*
 * TAlignment.cpp
 *
 *  Created on: 25.11.2011
 *      Author: bachmair
 */

#include "../include/TAlignment.hh"

TAlignment::TAlignment(int runNumber) {
	// TODO Auto-generated constructor stub
	sys = gSystem;
	stringstream  runString;
	runString.str("");
	runString<<runNumber;
	sys->MakeDirectory(runString.str().c_str());

	sys->cd(runString.str().c_str());
	stringstream  filepath;
	filepath.str("");
	filepath<<"clusterData."<<runNumber<<".root";
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
	initialiseHistos();
	cout<<"end initialise"<<endl;
}

TAlignment::~TAlignment() {
	// TODO Auto-generated destructor stub
	saveHistos();
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

void TAlignment::createVectors(UInt_t nEvents){
	int noHitDet=0;
	int falseClusterSizeDet=0;
	int noHitDia=0;
	int falseClusterSizeDia=0;
	int nCandidates=0;
	cout<<"ANALYSE VECTORS...."<<endl;
	for(UInt_t nEvent=0;nEvent<nEvents;nEvent++){
		TRawEventSaver::showStatusBar(nEvent,nEvents,100);
		eventReader->GetEvent(nEvent);
		//check if every plane has exactly one cluster
		bool candidate=true;
		for(UInt_t det=0;det<8&&candidate;det++){
			if(eventReader->getCluster()->at(det).size()!=1){
				candidate=false;
				noHitDet++;
				break;
			}
			if(candidate){
				TCluster cluster = eventReader->getCluster()->at(det).at(0);
				if(cluster.size()>2){
					candidate=false;
					falseClusterSizeDet++;
					break;
				}
			}
		}//for det
		if(candidate&&eventReader->getCluster()->at(8).size()!=1){
			candidate=false;
			cout<<"dia size:"<<eventReader->getCluster()->at(8).size()<<endl;
			noHitDia++;
		}
		if(candidate){
			TCluster cluster = eventReader->getCluster()->at(8).at(0);
			if(cluster.size()>2){
				candidate=false;
				falseClusterSizeDia++;
			}
		}
		//events with candidate=true areevents which have exactly one cluster in each plane
		// and the cluster size is 2
		if (candidate) nCandidates++;
	}
	cout<<"\n\nDetAnalysed "<<nEvents<<": "<<nCandidates<<" Candidates while "<< noHitDet<<" Events have not exactly one Cluster and "<<falseClusterSizeDet<<" Events have wrong cluster size"<<endl;
	cout<<"\n\nDiaAnalysed "<<nEvents-noHitDet-falseClusterSizeDet<<": "<<nCandidates<<" Candidates while "<< noHitDia<<" Events have not exactly one Cluster and "<<falseClusterSizeDia<<" Events have wrong cluster size"<<endl;
}

void TAlignment::initialiseHistos(){

}

void TAlignment::saveHistos(){

}
