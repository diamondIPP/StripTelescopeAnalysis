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
}

TClustering::~TClustering() {
	// TODO Auto-generated destructor stub
	delete eventReader;
	delete histSaver;
	sys->cd("..");
}

