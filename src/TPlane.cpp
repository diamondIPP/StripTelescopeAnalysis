//
//  TPlane.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TPlane.hh"
ClassImp(TPlane);
TPlane::TPlane(UInt_t planeNo,vector<TCluster> xClusters, vector<TCluster> yClusters,TPlaneProperties::enumDetectorType type) {
	verbosity=0;
	if(verbosity)cout<<"TPlane:"<<planeNo<<" xClusters:"<<xClusters.size()<<"\tyClusters:"<<yClusters.size()<<endl;
	for(UInt_t xCl=0;xCl<xClusters.size();xCl++)
		this->xClusters.push_back(xClusters.at(xCl));
	for(UInt_t yCl=0;yCl<yClusters.size();yCl++)
		this->yClusters.push_back(yClusters.at(yCl));
	this->type=type;
	this->planeNo=planeNo;

}

TPlane::TPlane(UInt_t planeNo,vector<TCluster> xClusters,TPlaneProperties::enumDetectorType type){
	verbosity=0;
	if(verbosity)cout<<"TPlane:"<<planeNo<<" xClusters:"<<xClusters.size()<<endl;
	for(UInt_t xCl=0;xCl<xClusters.size();xCl++)
		this->xClusters.push_back(xClusters.at(xCl));
	yClusters.clear();
	this->type=type;
	this->planeNo=planeNo;

}

TPlane::TPlane(const TPlane& rhs){
	this->verbosity=rhs.verbosity;
	if(verbosity>10)cout<<"Copy constructor of TPlane:"<<rhs.xClusters.size()<<" "<<rhs.yClusters.size()<<endl;
	for(UInt_t xCl=0;xCl<rhs.xClusters.size();xCl++){
		TCluster xCluster=rhs.xClusters.at(xCl);
		this->xClusters.push_back(xCluster);
	}
	for(UInt_t yCl=0;yCl<rhs.yClusters.size();yCl++){
		TCluster yCluster = rhs.yClusters.at(yCl);
		this->yClusters.push_back(yCluster);
	}
	this->type=rhs.type;
	this->planeNo=rhs.planeNo;
}

TPlane::~TPlane() {
	
}

enum TPlaneProperties::enumDetectorType TPlane::getDetectorType() const
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

void TPlane::setDetectorType(TPlaneProperties::enumDetectorType type)
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

/**
 * a Plane is valid if one and only one cluster is in each plane
 * @return
 */
bool TPlane::isValidPlane(){
	if(this->type == TPlaneProperties::kSilicon)
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

string TPlane::getDetectortypeString(TPlaneProperties::enumDetectorType type){
	switch (type){
	case TPlaneProperties::kSilicon: 	return "Silicon";
	case TPlaneProperties::kDiamond:	return "Diamond";
	default:		return "UNDEFINED";
	}
}

TCluster TPlane::getXCluster(UInt_t cl){
	if(cl<xClusters.size())
		return xClusters.at(cl);
	cerr<<"Xcluster does not exist....:"<<cl<<" "<<xClusters.size()<<endl;
	return TCluster();
}

TCluster TPlane::getYCluster(UInt_t cl){
	if(cl<yClusters.size())
		return yClusters.at(cl);
	cerr<<"Ycluster does not exist....:"<<cl<<" "<<yClusters.size()<<endl;
	return TCluster();
}

TCluster TPlane::getCluster(enumCoordinate cor, UInt_t cl){
	if(cor==X_COR)
		return getXCluster(cl);
	else
		return getYCluster(cl);
}

void TPlane::Print(UInt_t level)
{
	cout<< TCluster::Intent(level)<<getDetectortypeString(this->getDetectorType())<<"-Plane with "<<getNXClusters()<<"/"<<getNYClusters()<<endl;
	cout<<"X:";
	for(UInt_t i=0;i<getNXClusters();i++){
		this->xClusters.at(i).Print(level+1);
	}
	if(this->getDetectorType()==TPlaneProperties::kSilicon){
		cout<<"Y:";
		for(UInt_t i=0;i<getNYClusters();i++){
			this->yClusters.at(i).Print(level+1);
		}
	}
}
