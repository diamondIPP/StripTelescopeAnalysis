//
//  TEvent.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TEvent.hh"
ClassImp(TEvent);

TEvent::TEvent(){

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


}
