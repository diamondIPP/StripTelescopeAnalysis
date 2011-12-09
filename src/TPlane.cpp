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


