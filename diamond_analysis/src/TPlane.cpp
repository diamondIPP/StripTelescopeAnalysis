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
	this->verbosity=0;
	if(verbosity)cout<<"TPlane:"<<planeNo<<" xClusters:"<<xClusters.size()<<"\tyClusters:"<<yClusters.size()<<endl;
	this->xClusters.clear();
	this->yClusters.clear();
	for(UInt_t xCl=0;xCl<xClusters.size();xCl++)
		this->xClusters.push_back(xClusters.at(xCl));
	for(UInt_t yCl=0;yCl<yClusters.size();yCl++)
		this->yClusters.push_back(yClusters.at(yCl));
	this->type=type;
	this->planeNo=planeNo;

}

TPlane::TPlane(UInt_t planeNo,vector<TCluster> xClusters,TPlaneProperties::enumDetectorType type){
	this->verbosity=0;
	if(verbosity)cout<<"TPlane:"<<planeNo<<" xClusters:"<<xClusters.size()<<endl;
	this->xClusters.clear();
	for(UInt_t xCl=0;xCl<xClusters.size();xCl++)
		this->xClusters.push_back(xClusters.at(xCl));
	this->yClusters.clear();
	this->type=type;
	this->planeNo=planeNo;

}

/**
 * Copy Constructor for TPlane class
 * @param rhs
 */
TPlane::TPlane(const TPlane& rhs){
	this->verbosity=rhs.verbosity;
	this->xClusters.clear();
	this->yClusters.clear();
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

/**
 * Class Assignment function
 */
TPlane::TPlane& TPlane::operator =(const TPlane &src){

	type=src.type;
	planeNo=src.planeNo;
	verbosity=src.verbosity;
	xClusters.clear();
	for(UInt_t i=0;i<src.xClusters.size();i++)
		xClusters.push_back(src.xClusters.at(i));
	yClusters.clear();
	for(UInt_t i=0;i<src.yClusters.size();i++)
		yClusters.push_back(src.yClusters.at(i));
	return *this;
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


Float_t TPlane::getXPosition(UInt_t cl,TCluster::calculationMode_t mode,TH1F* histo){
	if(xClusters.size()>cl)
		return this->xClusters.at(cl).getPosition(mode,histo);
	else
		return N_INVALID;
}

Float_t TPlane::getYPosition(UInt_t cl,TCluster::calculationMode_t mode,TH1F* histo){
	if(yClusters.size()>cl)
		return this->yClusters.at(cl).getPosition(mode,histo);
	else
		return N_INVALID;
}

Float_t TPlane::getPosition(TPlaneProperties::enumCoordinate cor, UInt_t cl, TCluster::calculationMode_t mode,TH1F* histo)
{
	if(cor== TPlaneProperties::X_COR)
		return getXPosition(cl,mode,histo);
	else if(cor== TPlaneProperties::Y_COR)
		return getYPosition(cl,mode,histo);
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


TCluster TPlane::getXCluster(UInt_t cl){
	if(cl<xClusters.size())
		return ( xClusters.at(cl) );
	cerr<<"Xcluster does not exist....:"<<cl<<" "<<xClusters.size()<<endl;
	return ( TCluster() );
}

TCluster TPlane::getYCluster(UInt_t cl){
	if(cl<yClusters.size())
		return (yClusters.at(cl));
	cerr<<"Ycluster does not exist....:"<<cl<<" "<<yClusters.size()<<endl;
	return (TCluster());
}

TCluster TPlane::getCluster(TPlaneProperties::enumCoordinate cor, UInt_t cl){
	if(cor==TPlaneProperties::X_COR)
		return (getXCluster(cl));
	if(cor==TPlaneProperties::Y_COR)
		return (getYCluster(cl));
	else{
		cerr<<"Coordinate is neither X nor Y< return empty cluster"<<endl;
		return (TCluster());
	}
}

void TPlane::Print(UInt_t level)
{
	cout<< TCluster::Intent(level)<<TPlaneProperties::getDetectortypeString(this->getDetectorType())<<"-Plane with "<<getNXClusters()<<"/"<<getNYClusters()<<endl;
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
