/*
 * TTracking.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: bachmair
 */

#include "../include/TTracking.hh"

TTracking::TTracking(std::string pathName, std::string alignmentName):TADCEventReader(pathName){
	// TODO Auto-generated constructor stub
	alignmentFile=NULL;
	setAlignment(alignmentName);
	if(myAlignment!=NULL)
		myTrack=new TTrack(myAlignment);
	else
		myTrack=NULL;
}

TTracking::~TTracking() {
	// TODO Auto-generated destructor stub
}

bool TTracking::setAlignment(std::string alignmentName){
	if(this->alignmentFile!=NULL) alignmentFile->Delete();
	alignmentFile=NULL;
	cout<<"load AlignmentFile: \""<<alignmentName<<"\""<<endl;
	alignmentFile = new TFile(alignmentName.c_str());
	alignmentFile->GetObject("alignment",myAlignment);
	if(myAlignment==NULL){
		cerr<<"COULD NOT READTHE ALIGNMENT FILE..."<<endl;
		return false;
	}
	else{
		cout<<"Read the Alignment of AlignmentFile...\n"<<endl;
		myAlignment->PrintResults();
		return true;
	}
	return true;
}

bool TTracking::LoadEvent(UInt_t eventNumber){
	if(myTrack!=NULL){
		bool retVal=TADCEventReader::LoadEvent(eventNumber);
		if(retVal)
			myTrack->setEvent(this->getEvent());
		return retVal;
	}
	return false;
}
