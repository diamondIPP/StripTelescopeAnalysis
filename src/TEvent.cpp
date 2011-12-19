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
}

TEvent::~TEvent(){

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
		if(planes.at(plane).getDetectorType() == (TPlane::kSilicon)){
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

void TEvent::Print(UInt_t level){
	cout<<TCluster::Intent(level)<<"EventNo"<<getEventNumber()<<" with "<<getNPlanes()<< "Planes:"<<endl;
	for(UInt_t plane=0;plane<getNPlanes();plane++)
		planes.at(plane).Print(level+1);
}
