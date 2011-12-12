//
//  TPlane.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TPlane.hh"
ClassImp(TPlane);
TPlane::TPlane(vector<TCluster> xClusters, vector<TCluster> yClusters,enumDetectorType type) {
	this->xClusters = xClusters;
	this->yClusters = yClusters;
}

TPlane::~TPlane() {
	
}

UInt_t TPlane::getDetectorType() const
{
    return type;
}

void TPlane::setDetectorType(enumDetectorType type)
{
    this->type = type;
}


Float_t TPlane::getXPosition(UInt_t cl,TCluster::calculationMode_t mode){
	if(xClusters.size()>cl)
		return this->xClusters.at(cl).getPosition(mode);
	else
		return N_INVALID;
}

Float_t TPlane::getYPosition(UInt_t cl,TCluster::calculationMode_t mode){
	if(yClusters.size()>cl)
		return this->yClusters.at(cl).getPosition(mode);
	else
			return N_INVALID;
}


UInt_t TPlane::getNXClusters(){
	return xClusters.size();
}

UInt_t TPlane::getNYClusters(){
	return yClusters.size();
}

bool TPlane::isValidPlane(){
	if(this->type == kSilicon)
		return (getNXClusters()==1&&getNYClusters()==1);
	else
		return (getNXClusters()==1);
}
