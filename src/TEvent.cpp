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
}

TEvent::~TEvent(){

}
void TEvent::addPlane(TPlane plane, Int_t pos){
	if(pos==-1)
		this->planes.push_back(plane);
	else{
		if (planes.size()+1<pos){
			planes.resize(pos+1);
			planes.at(pos)=plane;
		}
		else
			planes.at(pos)=plane;
	}
}

bool TEvent::isValidSiliconTrack(){
	bool validTrack=true;
	for(UInt_t plane=0;plane<planes.size();plane++){
		if(planes.at(plane).getDetectorType() == (TPlane::kSilicon))
			validTrack=validTrack&&planes.at(plane).isValidPlane();
	}
	return validTrack;
}
