//
//  TPlane.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TPlane.hh"
ClassImp(TPlane);
TPlane::TPlane(vector<TCluster> xClusters, vector<TCluster> yClusters) {
	this->xClusters = xClusters;
	this->yClusters = yClusters;
}

TPlane::~TPlane() {
	
}
