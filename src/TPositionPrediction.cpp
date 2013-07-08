/*
 * TPositionPrediction.cpp
 *
 *  Created on: Dec 14, 2011
 *      Author: bachmair
 */

#include "../include/TPositionPrediction.hh"

TPositionPrediction::TPositionPrediction(Float_t zPos,Float_t zSigma,Float_t xPos,Float_t xSigma,Float_t xChi2, Float_t yPos, Float_t ySigma, Float_t yChi2,Float_t xPhi,Float_t yPhi) {
	// TODO Auto-generated constructor stub
	this->zPos=zPos;
	this->zSigma=zSigma;
	this->xPos=xPos;
	this->yPos=yPos;
	this->xSigma=xSigma;
	this->ySigma=ySigma;
	this->xChi2=xChi2;
	this->yChi2=yChi2;
	this->xPhi=xPhi;
	this->yPhi=yPhi;
	bValid=true;
}

TPositionPrediction::~TPositionPrediction() {
	// TODO Auto-generated destructor stub
}

void TPositionPrediction::setxChi2(Float_t chi2)
{
	xChi2 = chi2;
}

void TPositionPrediction::setxPos(Float_t pos)
{
	xPos = pos;
}

void TPositionPrediction::setxSigma(Float_t sigma)
{
	xSigma = sigma;
}


Float_t TPositionPrediction::getPosition(TPlaneProperties::enumCoordinate cor)
{
	switch (cor){
	case TPlaneProperties::X_COR: return this->getPositionX();break;
	case TPlaneProperties::Y_COR: return this->getPositionY();break;
	case TPlaneProperties::Z_COR: return this->getPositionZ();break;
	default: return N_INVALID;
	}
}

Float_t TPositionPrediction::getSigma(TPlaneProperties::enumCoordinate cor)
{
	switch (cor){
	case TPlaneProperties::X_COR: return this->getSigmaX();break;
	case TPlaneProperties::Y_COR: return this->getSigmaY();break;
	case TPlaneProperties::Z_COR: return this->getSigmaZ();break;
	default: return N_INVALID;
	}
}

/** gives the highest chi2 of both coordinates
 *
 * @return
 */
Float_t TPositionPrediction::getChi2(){
	return TMath::Max(yChi2,xChi2);
}

Float_t TPositionPrediction::getChi2(TPlaneProperties::enumCoordinate cor)
{
	switch (cor){
	case TPlaneProperties::X_COR: return this->getChi2X();break;
	case TPlaneProperties::Y_COR: return this->getChi2Y();break;
	default: return N_INVALID;
	}
}
Float_t TPositionPrediction::getPhi(TPlaneProperties::enumCoordinate cor){
	switch (cor){
	case TPlaneProperties::X_COR: return this->getPhiX();break;
	case TPlaneProperties::Y_COR: return this->getPhiY();break;
	default: return N_INVALID;
	}
}

void TPositionPrediction::setyChi2(Float_t chi2)
{
	yChi2 = chi2;
}

void TPositionPrediction::setyPos(Float_t pos)
{
	yPos = pos;
}

void TPositionPrediction::setySigma(Float_t sigma)
{
	ySigma = sigma;
}

void TPositionPrediction::Print(UInt_t level){
	cout<<TCluster::Intent(level)<<"TPositionPrediction::Print:"<<endl;
	cout<<TCluster::Intent(level+1)<<"Predicted Position of Track is: " <<xPos<<"+/-"<<xSigma<<" / "<<yPos<<"+/-"<<ySigma<<endl;
	cout<<TCluster::Intent(level+1)<< "The Fit has Chi2 values of "<<xChi2<<" and "<<yChi2<<endl;
}

