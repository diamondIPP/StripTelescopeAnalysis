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
	this->type=type;
}

TPlane::TPlane(vector<TCluster> xClusters,enumDetectorType type){
	this->xClusters = xClusters;
	this->type=type;
}

TPlane::~TPlane() {
	
}

enum TPlane::enumDetectorType TPlane::getDetectorType() const
{
	return type;
//	if(type==kSilicon)
//		return 1;
//	else if(type ==kDiamond)
//		return 2;
//	else if(type==kUndefined)
//		return 0;
//	else return 20;
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

Float_t TPlane::getPosition(enumCoordinate cor, UInt_t cl, TCluster::calculationMode_t mode)
{
	if(cor== X_COR)
		return getXPosition(cl,mode);
	else if(cor== Y_COR)
		return getYPosition(cl,mode);
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

string TPlane::getCoordinateString(enumCoordinate cor){
	switch (cor){
	case X_COR: return "X";break;
	case Y_COR: return "Y";break;
	case Z_COR: return "Z";break;
	case XY_COR:return "X&Y"; break;
	default: return "UNDEFINDED";
	}
}
