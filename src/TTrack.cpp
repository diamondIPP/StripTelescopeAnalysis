//
//  TTrack.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TTrack.hh"

TTrack::TTrack(TEvent event) {
	// TODO: check if event has the one & only one cluster flag
}

TTrack::~TTrack() {
	
}

void TTrack::setPositions(TEvent event) {
	for (int det = 0; det < 4; det++) {
		Float_t xOffset = this->getXOffset(det);
		Float_t yOffset = this->getYOffset(det);
		Float_t phiXOffset = this->getPhiXOffset(det);
		Float_t phiYOffset = this->getPhiXOffset(det);
		Float_t xPosition = event.getPlane(det).getXPosition(0);
		Float_t yPosition = event.getPlane(det).getYPosition(0);
		
		// apply offsets
		xPosition -= xOffset;
		yPosition -= yOffset;
		xPosition = (xPosition-128) * TMath::Cos(phiXOffset) - (yPosition-128) * TMath::Sin(-phiXOffset) + 128;
		yPosition = (xPosition-128) * TMath::Sin(-phiYOffset) + (yPosition-128) * TMath::Cos(phiYOffset) + 128;
		
		this->xPositions.push_back(xPosition);
		this->yPositions.push_back(yPosition);
	}
}