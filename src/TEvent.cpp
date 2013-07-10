//
//  TEvent.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TEvent.hh"
ClassImp(TEvent);

TEvent::TEvent(UInt_t nEvent){
	eventNumber=nEvent;
	planes.clear();
	verbosity=0;
}

TEvent::~TEvent(){

}

/**
 * COPY constructor
 * @param rhs right hand side of equal sign
 */
TEvent::TEvent(const TEvent& rhs){
	verbosity=rhs.verbosity;
	eventNumber = rhs.eventNumber;
	for(UInt_t pl=0;pl<rhs.planes.size();pl++)
		this->planes.push_back(rhs.planes.at(pl));
}

/**
 * class assignment function for TEvent Class
 * @param src
 * @return pointer to this TEvent
 */
TEvent &TEvent::operator=(const TEvent &src){
	planes.clear();
	for(UInt_t i=0;i<src.planes.size();i++)
		planes.push_back(src.planes.at(i));
	eventNumber=src.eventNumber;
	verbosity=src.verbosity;

	return (*this);
}

void TEvent::addPlane(TPlane plane, Int_t pos){
	if(pos==-1)
		this->planes.push_back(plane);
	else{
		if ((int)planes.size()<pos){
			planes.resize(pos+1);
			planes.at(pos)=plane;
		}
		else if((int)planes.size()==pos)
			planes.push_back(plane);
		else{
			planes.at(pos)=plane;
		}
	}
}

/**
 * checks if all silicon planes have one and only one cluster in each detector layer
 */
bool TEvent::isValidSiliconEvent(){
	bool validTrack=true;
	for(UInt_t plane=0;plane<planes.size();plane++){
		if(planes.at(plane).getDetectorType() == (TPlaneProperties::kSilicon)){
			validTrack=validTrack&&planes.at(plane).isValidPlane();
		}
	}
	return validTrack;
}

void TEvent::setVerbosity(UInt_t verbosity)
{
	this->verbosity=verbosity;
}

bool TEvent::isMasked(){
	//todo
	return false;
}

UInt_t TEvent::getNPlanes(){
	return planes.size();
}
UInt_t TEvent::getNClusters(UInt_t det){
	TPlaneProperties::enumCoordinate cor;
	if(det%2==0)
		cor=TPlaneProperties::X_COR;
	else
		cor=TPlaneProperties::Y_COR;
	UInt_t plane=det/2;
	if(verbosity>12)cout<<"TEvent::getNclusters of det "<<det<<" <=> plane: "<<plane<<" "<<TPlaneProperties::getCoordinateString(cor)<<endl;
	if(det%2==0)
		return getNXClusters(plane);
	else return getNYClusters(plane);

}
UInt_t TEvent::getNXClusters(UInt_t plane){
	if(verbosity>10)
		cout<<"TEvent::getNXClusters of plane "<<plane<<endl;
	if(plane<planes.size())
		return planes.at(plane).getNXClusters();
	else
		return 0;
}

UInt_t TEvent::getNYClusters(UInt_t plane){
	if(verbosity>2)
		cout<<"TEvent::getNYClusters of plane "<<plane<<endl;
	if(plane<planes.size())
		return planes.at(plane).getNYClusters();
	else
		return 0;
}

TCluster TEvent::getCluster(UInt_t plane,TPlaneProperties::enumCoordinate cor, UInt_t cl){
	if (plane<planes.size())
		return planes.at(plane).getCluster(cor,cl);
	cerr<< "Plane does not exist: "<<plane<<" "<<planes.size()<<endl;
	return TCluster();
}

TCluster TEvent::getCluster(UInt_t det,UInt_t cl){
	UInt_t plane = det/2;
	TPlaneProperties::enumCoordinate cor;
	if (det%2==0) cor =TPlaneProperties::X_COR;
	else cor = TPlaneProperties::Y_COR;
	return this->getCluster(plane,cor,cl);
}

UInt_t TEvent::getClusterSize(UInt_t det, UInt_t cl){
	TCluster cluster = getCluster(det,cl);
	if(verbosity>20)cluster.Print();
	return cluster.size();
}

UInt_t TEvent::getClusterSeedSize(UInt_t det, UInt_t cl){
	TCluster cluster = getCluster(det,cl);
	if(verbosity>20)cluster.Print();
	return cluster.seedSize();
}

UInt_t TEvent::getClusterSize(UInt_t plane,TPlaneProperties::enumCoordinate cor,UInt_t cl){
	return getCluster(plane,cor,cl).size();
}

Float_t TEvent::getPosition(UInt_t det, UInt_t cl,TCluster::calculationMode_t mode,TH1F *histo)
{
	TCluster cluster = getCluster(det,cl);
	return cluster.getPosition(mode,histo);
}

bool TEvent::hasInvalidReadout(){
	bool invalidReadout=false;
//	cout<<eventNumber<<" TEvent: invalidReadout"<<endl;
	for(UInt_t plane=0;plane<getNPlanes()&&!invalidReadout;plane++){
		invalidReadout = planes.at(plane).hasInvalidReadout()||invalidReadout;
	}
	return invalidReadout;
}
void TEvent::Print(UInt_t level){
	cout<<TCluster::Intent(level)<<"EventNo"<<getEventNumber()<<" with "<<getNPlanes()<< "Planes:"<<endl;
	for(UInt_t plane=0;plane<getNPlanes();plane++)
		planes.at(plane).Print(level+1);
}
