/*
 * TResidual.cpp
 *
 *  Created on: Jan 11, 2012
 *      Author: bachmair
 */

#include "../include/TResidual.hh"

using namespace std;
TResidual::TResidual() {
	resXMean=0;
	resXSigma=0;
	resYMean=0;
	resYSigma=0;
	sumRx=0;
	sumRy=0;
	sumVx=0;
	sumVy=0;
	sumV2x=0;
	sumV2y=0;
	sumVRx=0;
	sumVRy=0;
	nUsedTracks=0;
	bTestResidual=false;
}

TResidual::~TResidual() {
	// TODO Auto-generated destructor stub
}

void TResidual::Print(UInt_t level){
	cout<<"\n\n";
	cout<<TCluster::Intent(level)<<"residual results of "<<nUsedTracks<<" used Tracks"<<endl;
	cout<<TCluster::Intent(level)<<"\tX: "<<getXMean()<<"+/-"<<getXSigma()<<"\t"<<sumRx<<" "<<sumVx<<" "<<sumV2x<<" "<<sumVRx<<endl;
	cout<<TCluster::Intent(level)<<"\tY: "<<getYMean()<<"+/-"<<getYSigma()<<"\t"<<sumRy<<" "<<sumVy<<" "<<sumV2y<<" "<<sumVRy<<endl;
	cout<<TCluster::Intent(level)<<"\t Xoff: "<<getXOffset()<<"\tPhiXoff: "<<getPhiXOffset()<<endl;
	cout<<TCluster::Intent(level)<<"\t Yoff: "<<getYOffset()<<"\tPhiYoff: "<<getPhiYOffset()<<endl;
	cout<<"\n\n\n"<<endl;
}

void TResidual::addDataPoint(Float_t deltaX, Float_t predX, Float_t deltaY, Float_t predY)
{
	this->resYMean += deltaY;
	this->resYSigma += deltaY*deltaY;
	this->resXMean += deltaX;
	this->resXSigma +=deltaX*deltaX;
	this->sumRx+=deltaX;
	this->sumRy+=deltaY;
	this->sumVx+=predY;
	this->sumVy+=predX;
	this->sumV2x+=predY*predY;//todo check
	this->sumV2y+=predX*predX;//todo check
	this->sumVRx+=predY*deltaX;
	this->sumVRy+=predX*deltaY;
	this->nUsedTracks++;
}

Float_t TResidual::getXSigma()
{
	if(bTestResidual) return 100000;
	if (nUsedTracks!=0)
		return TMath::Sqrt(this->resXSigma / (Double_t)this->nUsedTracks - getXMean()*getXMean());
	return N_INVALID;
}



Float_t TResidual::getYSigma()
{
	if(bTestResidual) return 100000;
	if (nUsedTracks!=0)
		return TMath::Sqrt(this->resYSigma / (Double_t)this->nUsedTracks - getYMean()*getYMean());
	return N_INVALID;
}

Float_t TResidual::getXOffset()
{
	Float_t variableDif= (nUsedTracks * sumV2x - sumVx * sumVx);
	if (variableDif!=0)
		return (sumRx * sumV2x - sumVRx * sumVx) / variableDif;
	return N_INVALID;
}



Float_t TResidual::getYOffset()
{
	Float_t variableDif= (nUsedTracks * sumV2y - sumVy * sumVy);
	if (variableDif!=0)
		return (sumRy * sumV2y - sumVRy * sumVy) / variableDif;
	return N_INVALID;
}

Float_t TResidual::getPhiXOffset()
{
	Float_t variableDif = (nUsedTracks * sumV2x - sumVx * sumVx);
	if(variableDif!=0)
		return TMath::ATan(-(nUsedTracks * sumVRx - sumRx * sumVx) / variableDif);
	return N_INVALID;
}



Float_t TResidual::getXMean()
{
	if(bTestResidual)
		return 0;
	if(nUsedTracks!=0)
		return  (this->resXMean/(Double_t)this->nUsedTracks);
	return N_INVALID;
}

Float_t TResidual::getYMean()
{
	if(bTestResidual)
		return 0;
	if(nUsedTracks!=0)
		return  (this->resYMean/(Double_t)this->nUsedTracks);
	return N_INVALID;
}

Float_t TResidual::getPhiYOffset()
{
	Float_t variableDif=(nUsedTracks * sumV2y - sumVy * sumVy);
	if(variableDif!=0)
		return TMath::ATan((nUsedTracks* sumVRy - sumRy * sumVy) / variableDif);
	return N_INVALID;
}








