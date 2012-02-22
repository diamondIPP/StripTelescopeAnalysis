/*
 * TTracking.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: bachmair
 */

#include "../include/TTracking.hh"

TTracking::TTracking(std::string pathName, std::string alignmentName,UInt_t runNumber):TADCEventReader(pathName,runNumber){
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

TPositionPrediction *TTracking::predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPrint)
{
	if(myTrack==0)
		return 0;
	return myTrack->predictPosition(subjectPlane,vecRefPlanes,bPrint);
}

Float_t TTracking::getXPosition(UInt_t plane)
{
	if(myTrack==0)
			return 0;
		return myTrack->getXPosition(plane);
}

Float_t TTracking::getYPosition(UInt_t plane)
{
	if(myTrack==0)
			return 0;
		return myTrack->getYPosition(plane);
}

Float_t TTracking::getZPosition(UInt_t plane)
{
	if(myTrack==0)
			return 0;
		return myTrack->getZPosition(plane);
}

Float_t TTracking::getMeasured(TPlane::enumCoordinate cor, UInt_t plane,TCluster::calculationMode_t mode)
{
	if (myTrack == 0)
		return 0;
	return myTrack->getMeasured(cor,plane,mode);
}

bool TTracking::setAlignment(std::string alignmentName){
	if(this->alignmentFile!=NULL) alignmentFile->Delete();
	alignmentFile=NULL;
	cout<<"load AlignmentFile: \""<<alignmentName<<"\""<<endl;
	alignmentFile = new TFile(alignmentName.c_str());
	alignmentFile->GetObject("alignment",myAlignment);
	if(myAlignment==NULL){
		cerr<<"COULD NOT READ THE ALIGNMENT FILE..."<<endl;
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

Float_t  TTracking::getStripXPositionOfCluster(UInt_t plane,TCluster xCluster, Float_t yPred,TCluster::calculationMode_t mode){
	if(myTrack==0)
				return 0;
			return myTrack->getStripXPositionOfCluster(plane,xCluster,yPred,mode);
}
Float_t  TTracking::getStripXPosition(UInt_t plane,Float_t yPred,TCluster::calculationMode_t mode){
	if(myTrack==0)
				return 0;
			return myTrack->getStripXPosition(plane,yPred,mode);

}
Float_t  TTracking::getPositionOfCluster(TPlane::enumCoordinate cor,UInt_t plane,TCluster xCluster,TCluster yCluster, TCluster::calculationMode_t mode){
	if(myTrack==0)
				return 0;
			return myTrack->getPositionOfCluster(cor,plane,xCluster,yCluster,mode);
}
Float_t TTracking::getPositionOfCluster(UInt_t det, TCluster cluster, Float_t predictedPerpPosition, TCluster::calculationMode_t mode){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionOfCluster(det,cluster,predictedPerpPosition,mode);
}
Float_t  TTracking::getPosition(TPlane::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode){
	if(myTrack==0)
		return 0;
	return myTrack->getPosition(cor,plane,mode);
}
Float_t TTracking::getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionInDetSystem(det,xPred,yPred);
}
