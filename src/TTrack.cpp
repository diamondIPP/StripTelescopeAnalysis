//
//  TTrack.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TTrack.hh"

TTrack::TTrack(TDetectorAlignment alignment) {
	// TODO: check if event has the one & only one cluster flag
	verbosity=0;
	this->alignment = alignment;
	event=NULL;
}

TTrack::~TTrack() {
	
}

//void TTrack::setPositions(TEvent event) {
//	for (int plane = 0; plane < 4; plane++) {
//		Float_t xOffset = this->getXOffset(plane);
//		Float_t yOffset = this->getYOffset(plane);
//		Float_t phiXOffset = this->getPhiXOffset(plane);
//		Float_t phiYOffset = this->getPhiXOffset(plane);
//		Float_t xPosition = event.getPlane(plane).getXPosition(0);
//		Float_t yPosition = event.getPlane(plane).getYPosition(0);
//		
//		// apply offsets
//		xPosition -= xOffset;
//		yPosition -= yOffset;
//		xPosition = (xPosition-128) * TMath::Cos(phiXOffset) - (yPosition-128) * TMath::Sin(-phiXOffset) + 128;
//		yPosition = (xPosition-128) * TMath::Sin(-phiYOffset) + (yPosition-128) * TMath::Cos(phiYOffset) + 128;
//		
//		this->xPositions.push_back(xPosition);
//		this->yPositions.push_back(yPosition);
//	}
//}

UInt_t TTrack::getNClusters(int det) {
	if(event==NULL)return N_INVALID;
	if (det%2 == 0) {
		int plane = det / 2;
		return this->event->getPlane(plane).getNXClusters();
	}
	else {
		int plane = (det-1) / 2;
		return this->event->getPlane(plane).getNYClusters();
	}
}

Float_t TTrack::getXPosition(UInt_t plane) {
	if(event==NULL)return N_INVALID;
	if(event->getNXClusters(plane)!=1||event->getNYClusters(plane)!=1)
		return N_INVALID;
	// get offsets
	Float_t xOffset = this->getXOffset(plane);
	Float_t yOffset = this->getYOffset(plane);
	Float_t phiXOffset = this->getPhiXOffset(plane);
//	Float_t phiYOffset = this->getPhiXOffset(plane);
	Float_t xPosition = event->getPlane(plane).getXPosition(0);
	Float_t yPosition = event->getPlane(plane).getYPosition(0);
	
	// apply offsets
	xPosition -= xOffset;
	yPosition -= yOffset;
	xPosition = (xPosition-128) * TMath::Cos(phiXOffset) - (yPosition-128) * TMath::Sin(-phiXOffset) + 128;
//	yPosition = (xPosition-128) * TMath::Sin(-phiYOffset) + (yPosition-128) * TMath::Cos(phiYOffset) + 128;
	
	return xPosition;
}

Float_t TTrack::getYPosition(UInt_t plane) {
	if(event==NULL)return N_INVALID;
	if(event->getNXClusters(plane)!=1||event->getNYClusters(plane)!=1)
		return N_INVALID;
	// get offsets
	Float_t xOffset = this->getXOffset(plane);
	Float_t yOffset = this->getYOffset(plane);
	Float_t phiXOffset = this->getPhiXOffset(plane);
	Float_t phiYOffset = this->getPhiXOffset(plane);
	Float_t xPosition = event->getPlane(plane).getXPosition(0);
	Float_t yPosition = event->getPlane(plane).getYPosition(0);
	
	// apply offsets
	xPosition -= xOffset;
	yPosition -= yOffset;
	xPosition = (xPosition-128) * TMath::Cos(phiXOffset) - (yPosition-128) * TMath::Sin(-phiXOffset) + 128;
	yPosition = (xPosition-128) * TMath::Sin(-phiYOffset) + (yPosition-128) * TMath::Cos(phiYOffset) + 128;
	
	return yPosition;
}

Float_t TTrack::getPosition(TPlane::enumCoordinate cor, UInt_t plane){
	if(event==NULL)return 0;
	if(cor == TPlane::X_COR)
		return getXPosition(plane);
	else if(cor == TPlane::Y_COR)
		return getYPosition(plane);
	else
		return N_INVALID;
}


TPositionPrediction TTrack::predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes)
{
	if(event==NULL){
		cerr<<"TTrack:predictPosition no ReferencePlanes are defined..."<<endl;
		TPositionPrediction prediction;
		return prediction;
	}
	if(vecRefPlanes.size()==0){
		cerr<<"TTrack:predictPosition no ReferencePlanes are defined..."<<endl;
		TPositionPrediction prediction;
		return prediction;
	}
	if(vecRefPlanes.size()==1){
		if(verbosity>3)	cout<<"TTrack::predictPosition with 1 refPlane"<<endl;
		TPositionPrediction prediction(getXPosition(vecRefPlanes.at(0)), 0.,0.,getYPosition(vecRefPlanes.at(0)),0.,0.);
		return prediction;
	}
	TLinearFitter* linFitX = new TLinearFitter(1,"pol1","D");
new TLinearFitter((Int_t)1,"pol1","D");
	TLinearFitter* linFitY= new TLinearFitter(1,"pol1");
	vector<Double_t>vecXPos;
	vector<Double_t>vecYPos;
	vector<Double_t>vecZPos;
	vector<Double_t> vecZSigma;
	vector<Double_t> zPosVec;
	for(UInt_t pl=0;pl<vecRefPlanes.size();pl++){
		UInt_t plane=vecRefPlanes.at(pl);
		vecXPos.push_back(getXPosition(plane));
		vecYPos.push_back(getYPosition(plane));
		vecZPos.push_back(alignment.GetZOffset(plane));
		vecZSigma.push_back(0.0);
		zPosVec.clear();
		zPosVec.push_back(alignment.GetZOffset(plane));
		linFitX->AddPoint(&zPosVec.at(0),(Double_t)getXPosition(plane),0.001);
		linFitY->AddPoint(&zPosVec.at(0),(Double_t)getYPosition(plane),0.001);
		if(verbosity>3)	cout<<"Add "<<getXPosition(plane)<<"/"<<getYPosition(plane)<<"/"<<alignment.GetZOffset(plane)<<endl;
	}
//	linFitX->AssignData(vecXPos.size(),1,&vecZPos.at(0),&vecXPos.at(0),&vecZPos.at(0));
//	linFitY->AssignData(vecYPos.size(),1,&vecZPos.at(0),&vecYPos.at(0),&vecZPos.at(0));
	linFitX->Eval();
	linFitY->Eval();
	Float_t zPos = alignment.GetZOffset(subjectPlane);
	Float_t zSigma = 0; //todo
	Float_t mx = linFitX->GetParameter(1);
	Float_t sigma_mx = linFitX->GetParError(1);
	Float_t bx = linFitX->GetParameter(0);
	Float_t sigma_bx = linFitX->GetParError(0);
	Float_t my = linFitY->GetParameter(1);
	Float_t sigma_my = linFitY->GetParError(1);
	Float_t by = linFitY->GetParameter(0);
	Float_t sigma_by = linFitY->GetParError(0);
	Float_t xChi2 = linFitX->GetChisquare();
	Float_t yChi2 = linFitY->GetChisquare();
	Float_t xPos = mx*zPos+bx;
	Float_t yPos = my*zPos+by;
	Float_t xSigma = (zPos*sigma_mx)*(zPos*sigma_mx)+(mx*zSigma)*(mx*zSigma)+(sigma_bx*sigma_bx);
	xSigma = TMath::Sqrt(xSigma);
	Float_t ySigma = (zPos*sigma_my)*(zPos*sigma_my)+(my*zSigma)*(my*zSigma)+(sigma_by*sigma_by);
	ySigma = TMath::Sqrt(ySigma);
	TPositionPrediction prediction(xPos,xSigma,xChi2,yPos,ySigma,yChi2);
	if(verbosity>3)	cout<<mx<<"+/-"<<sigma_mx<<"    "<<bx<<"+/-"<<sigma_bx<<endl;
	if(verbosity>3)	cout<<my<<"+/-"<<sigma_my<<"    "<<by<<"+/-"<<sigma_by<<endl;
	return prediction;
}

vector<Float_t> TTrack::getSiXPositions() {
	vector<Float_t> xPositions;
	if(event==NULL)return xPositions;
	for (int plane = 0; plane < 4; plane++) {
		xPositions.push_back(this->getXPosition(plane));
	}
	return xPositions;
}

vector<Float_t> TTrack::getSiYPositions() {

	vector<Float_t> yPositions;
	if(event==NULL)return yPositions;
	for (int plane; plane < 4; plane++) {
		yPositions.push_back(this->getYPosition(plane));
	}
	return yPositions;
}

void TTrack::setDetectorAlignment(TDetectorAlignment alignment)
{
	this->alignment=alignment;
}


