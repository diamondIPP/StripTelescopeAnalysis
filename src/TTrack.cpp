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
	
	this->alignment = alignment;
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
	if (det%2 == 0) {
		int plane = det / 2;
		return this->event->getPlane(plane).getNXClusters();
	}
	else {
		int plane = (det-1) / 2;
		return this->event->getPlane(plane).getNYClusters();
	}
}

Float_t TTrack::getXPosition(int plane) {
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

Float_t TTrack::getYPosition(int plane) {
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

vector<Float_t> TTrack::getSiXPositions() {
	vector<Float_t> xPositions;
	for (int plane = 0; plane < 4; plane++) {
		xPositions.push_back(this->getXPosition(plane));
	}
	return xPositions;
}

vector<Float_t> TTrack::getSiYPositions() {
	vector<Float_t> yPositions;
	for (int plane; plane < 4; plane++) {
		yPositions.push_back(this->getYPosition(plane));
	}
	return yPositions;
}