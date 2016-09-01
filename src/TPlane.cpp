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
	this->SetClusters(xClusters, yClusters);
	this->type=type;
	this->planeNo=planeNo;
//	this->adc.clear();
//	this->ped.clear();
//	this->pedCMN.clear();
//	this->pedSigma.clear();
//	this->pedSigmaCMN.clear();
//	this->rawSignal.clear();
//	this->rawSignalCMN.clear();
//	this->cmNoise = 0;
}

TPlane::TPlane(UInt_t planeNo,vector<TCluster> xClusters,TPlaneProperties::enumDetectorType type) {
	this->verbosity = 0;
	if (verbosity)cout << "TPlane:" << planeNo << " xClusters:" << xClusters.size() << endl;
	this->SetXClusters(xClusters);
	this->yClusters.clear();
	this->type = type;
	this->planeNo = planeNo;
//	this->adc.clear();
//	this->ped.clear();
//	this->pedCMN.clear();
//	this->pedSigma.clear();
//	this->pedSigmaCMN.clear();
//	this->rawSignal.clear();
//	this->rawSignalCMN.clear();
//	this->cmNoise = 0;
}

//TPlane::TPlane(UInt_t planeNo,vector<TCluster> xClusters,TPlaneProperties::enumDetectorType type, Float_t ecmNoise){
//	this->verbosity=0;
//	if(verbosity)cout<<"TPlane:"<<planeNo<<" xClusters:"<<xClusters.size()<<endl;
//	this->SetXClusters(xClusters);
//	this->yClusters.clear();
//	this->type=type;
//	this->planeNo=planeNo;
//	this->adc.clear();
//	this->ped.clear();
//	this->pedCMN.clear();
//	this->pedSigma.clear();
//	this->pedSigmaCMN.clear();
//	this->rawSigna.clear();
//	this->rawSignalCMN.clear();
//	this->cmNoise = ecmNoise;
//}

//TPlane::TPlane(UInt_t planeNo,vector<TCluster> xCluster, TPlaneProperties::enumDetectorType type, UShort_t *eadc, Float_t *eped, Float_t *epedCMN, Float_t *epedSigma, Float_t *epedSigmaCMN, Float_t *erawSignal, Float_t *erawSignalCMN, Float_t ecmNoise){
//	this->verbosity=0;
//	if(verbosity)cout<<"TPlane:"<<planeNo<<" xClusters:"<<xClusters.size()<<endl;
//	this->SetXClusters(xClusters);
//	this->yClusters.clear();
//	this->type=type;
//	this->planeNo=planeNo;
//	this->SetSignalValues(eadc, eped, epedCMN, epedSigma, epedSigmaCMN, erawSignal, erawSignalCMN, ecmNoise);
//}

//TPlane::TPlane(UInt_t planeNo,vector<TCluster> xClusters,TPlaneProperties::enumDetectorType type)

/**
 * Copy Constructor for TPlane class
 * @param rhs
 */
TPlane::TPlane(const TPlane& rhs){

	this->verbosity=rhs.verbosity;
	this->xClusters.clear();
	this->yClusters.clear();
//	this->adc.clear();
//	this->ped.clear();
//	this->pedCMN.clear();
//	this->pedSigma.clear();
//	this->pedSigmaCMN.clear();
//	this->rawSignal.clear();
//	this->rawSignalCMN.clear();
//	this->cmNoise = 0;
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
//	for(UInt_t iadc = 0; iadc<rhs.adc.size(); iadc++) {
//		this->adc.push_back(rhs.adc.at(iadc));
//	}
//	for(UInt_t iped = 0; iped<rhs.ped.size(); iped++) {
//		this->ped.push_back(rhs.ped.at(iped));
//	}
//	for(UInt_t ipedCMN = 0; ipedCMN<rhs.pedCMN.size(); ipedCMN++) {
//		this->pedCMN.push_back(rhs.pedCMN.at(ipedCMN));
//	}
//	for(UInt_t ipedSigma = 0; ipedSigma<rhs.pedSigma.size(); ipedSigma++) {
//		this->pedSigma.push_back(rhs.pedSigma.at(ipedSigma));
//	}
//	for(UInt_t ipedSigmaCMN = 0; ipedSigmaCMN<rhs.pedSigmaCMN.size(); ipedSigmaCMN++) {
//		this->pedSigmaCMN.push_back(rhs.pedSigmaCMN.at(ipedSigmaCMN));
//	}
//	for(UInt_t irawSignal = 0; irawSignal<rhs.rawSignal.size(); irawSignal++) {
//		this->rawSignal.push_back(rhs.rawSignal.at(irawSignal));
//	}
//	for(UInt_t irawSignalCMN = 0; irawSignalCMN<rhs.rawSignalCMN.size(); irawSignalCMN++) {
//		this->rawSignalCMN.push_back(rhs.rawSignalCMN.at(irawSignalCMN));
//	}
//	this->cmNoise = rhs.cmNoise;
}

/**
 * Class Assignment function
 */
TPlane& TPlane::operator =(const TPlane &src){
    if (this == &src)
        return *this;
    TPlane::operator=(src);
	type=src.type;
	planeNo=src.planeNo;
	verbosity=src.verbosity;
	xClusters.clear();
//	adc.clear();
//	ped.clear();
//	pedCMN.clear();
//	pedSigma.clear();
//	pedSigmaCMN.clear();
//	rawSignal.clear();
//	rawSignalCMN.clear();
	for(UInt_t i=0;i<src.xClusters.size();i++)
		xClusters.push_back(src.xClusters.at(i));
	yClusters.clear();
	for(UInt_t i=0;i<src.yClusters.size();i++)
		yClusters.push_back(src.yClusters.at(i));
//	for(UInt_t iadc = 0; iadc<src.adc.size(); iadc++) {
//		adc.push_back(src.adc.at(iadc));
//	}
//	for(UInt_t iped = 0; iped<src.ped.size(); iped++) {
//		ped.push_back(src.ped.at(iped));
//	}
//	for(UInt_t ipedCMN = 0; ipedCMN<src.pedCMN.size(); ipedCMN++) {
//		pedCMN.push_back(src.pedCMN.at(ipedCMN));
//	}
//	for(UInt_t ipedSigma = 0; ipedSigma<src.pedSigma.size(); ipedSigma++) {
//		pedSigma.push_back(src.pedSigma.at(ipedSigma));
//	}
//	for(UInt_t ipedSigmaCMN = 0; ipedSigmaCMN<src.pedSigmaCMN.size(); ipedSigmaCMN++) {
//		pedSigmaCMN.push_back(src.pedSigmaCMN.at(ipedSigmaCMN));
//	}
//	for(UInt_t irawSignal = 0; irawSignal<src.rawSignal.size(); irawSignal++) {
//		rawSignal.push_back(src.rawSignal.at(irawSignal));
//	}
//	for(UInt_t irawSignalCMN = 0; irawSignalCMN<src.rawSignalCMN.size(); irawSignalCMN++) {
//		rawSignalCMN.push_back(src.rawSignalCMN.at(irawSignalCMN));
//	}
//	cmNoise = src.cmNoise;
	return *this;
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


Float_t TPlane::getXPosition(UInt_t cl,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo){
	if(xClusters.size()>cl)
		return this->xClusters.at(cl).getPosition(cmnCorrected,mode,histo);
	else
		return N_INVALID;
}

Float_t TPlane::getYPosition(UInt_t cl,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo){
	if(yClusters.size()>cl)
		return this->yClusters.at(cl).getPosition(cmnCorrected,mode,histo);
	else
		return N_INVALID;
}

Float_t TPlane::getPosition(TPlaneProperties::enumCoordinate cor, UInt_t cl, bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo)
{
	if(cor== TPlaneProperties::X_COR)
		return getXPosition(cl,cmnCorrected,mode,histo);
	else if(cor== TPlaneProperties::Y_COR)
		return getYPosition(cl,cmnCorrected,mode,histo);
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
		cerr<<"Plane "<<planeNo<<": Coordinate is neither X nor Y< return empty cluster "<<cor<<endl;
		return (TCluster());
	}
}


void TPlane::SetClusters(vector<TCluster> xClusters,
		vector<TCluster> yClusters) {
	this->SetXClusters(xClusters);
	this->SetYClusters(yClusters);
}

void TPlane::SetXClusters(vector<TCluster> xClusters) {
	this->xClusters.clear();
	for(UInt_t xCl=0;xCl<xClusters.size();xCl++)
		this->xClusters.push_back(xClusters.at(xCl));
}

void TPlane::SetYClusters(vector<TCluster> yClusters) {
	this->yClusters.clear();
	for(UInt_t yCl=0;yCl<yClusters.size();yCl++)
		this->yClusters.push_back(yClusters.at(yCl));

}


bool TPlane::hasInvalidReadout(){
//	cout<<"\tTPlane::hasInvalidReadout "<<planeNo<<": "<<xClusters.size()<<"-"<<yClusters.size()<<endl;
	bool invalidReadout=false;
//	cout<<"\tX:\t"<<endl;
	UInt_t xCl;
	for(xCl=0;xCl<xClusters.size()&&!invalidReadout;xCl++)
		invalidReadout = invalidReadout || xClusters.at(xCl).hasInvalidReadout();
//	if(invalidReadout)
//		cout<<planeNo<<" X"<<xCl<<endl;
//	cout<<"\tY:"<<endl;
	UInt_t yCl;
	for(yCl=0;yCl<yClusters.size()&&!invalidReadout;yCl++){
		invalidReadout = invalidReadout || yClusters.at(yCl).hasInvalidReadout();
//		if(invalidReadout)
//			cout<<planeNo<<"Y"<<xCl<<endl;
	}//	cout<<"\treturning: "<<invalidReadout<<endl;
	return invalidReadout;
}

//void TPlane::SetSignalValues(UShort_t *eadc, Float_t *eped, Float_t *epedCMN, Float_t *epedSigma, Float_t *epedSigmaCMN, Float_t *erawSignal, Float_t *erawSignalCMN, Float_t ecmNoise){
//	this->adc.clear();
//	this->ped.clear();
//	this->pedCMN.clear();
//	this->pedSigma.clear();
//	this->pedSigmaCMN.clear();
//	this->rawSignal.clear();
//	this->rawSignalCMN.clear();
//	for (UInt_t ch = 0; ch < TPlaneProperties::getNChannelsDiamond(); ch++){
//		this->adc.push_back(eadc[ch]);
//		this->ped.push_back(eped[ch]);
//		this->pedCMN.push_back(epedCMN[ch]);
//		this->pedSigma.push_back(epedSigma[ch]);
//		this->pedSigmaCMN.push_back(epedSigmaCMN[ch]);
//		this->rawSignal.push_back(erawSignal[ch]);
//		this->rawSignalCMN.push_back(erawSignalCMN[ch]);
//	}
//	this->cmNoise = ecmNoise;
//}

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
