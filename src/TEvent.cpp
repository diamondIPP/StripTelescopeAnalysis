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

TEvent::TEvent(const TEvent& rhs){
	verbosity=rhs.verbosity;
	eventNumber = rhs.eventNumber;
	for(UInt_t pl=0;pl<rhs.planes.size();pl++)
		this->planes.push_back(rhs.planes.at(pl));
}

void TEvent::addPlane(TPlane plane, Int_t pos){
	if(pos==-1)
		this->planes.push_back(plane);
	else{
		if (planes.size()<pos){
			planes.resize(pos+1);
			planes.at(pos)=plane;
		}
		else if(planes.size()==pos)
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

bool TEvent::isMasked(){
	//todo
	return false;
}

UInt_t TEvent::getNPlanes(){
	return planes.size();
}
UInt_t TEvent::getNClusters(UInt_t det){
	TPlane::enumCoordinate cor;
	UInt_t plane=det/2;
	if(det%2==0)
		return getNXClusters(plane);
	else return getNYClusters(plane);

}

TCluster TEvent::getCluster(UInt_t plane,TPlane::enumCoordinate cor, UInt_t cl){
	if (plane<planes.size())
		return planes.at(plane).getCluster(cor,cl);
	else
		return TCluster();
}

TCluster TEvent::getCluster(UInt_t det,UInt_t cl){
	UInt_t plane = det/2;
	TPlane::enumCoordinate cor;
	if (det%2==0) cor =TPlane::X_COR;
	else cor = TPlane::Y_COR;
	return this->getCluster(plane,cor,cl);
}

UInt_t TEvent::getClusterSize(UInt_t det, UInt_t cl){
	TCluster cluster = getCluster(det,cl);
	if(verbosity>20)cluster.Print();
	return cluster.size();
}

UInt_t TEvent::getClusterSize(UInt_t plane,TPlane::enumCoordinate cor,UInt_t cl){
	return getCluster(plane,cor,cl).size();
}

void TEvent::Print(UInt_t level){
	cout<<TCluster::Intent(level)<<"EventNo"<<getEventNumber()<<" with "<<getNPlanes()<< "Planes:"<<endl;
	for(UInt_t plane=0;plane<getNPlanes();plane++)
		planes.at(plane).Print(level+1);
}
