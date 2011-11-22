/*
 * TCluster.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#include "../include/TCluster.hh"

TCluster::TCluster(int eventNumber,int seedSigma,int hitSigma) {
	// TODO Auto-generated constructor stub
	numberOfSeeds=0;
	numberOfHits=0;
	this->seedSigma=seedSigma;
	this->hitSigma=hitSigma;
}

TCluster::~TCluster() {
	// TODO Auto-generated destructor stub
}

void TCluster::addChannel(int ch, Float_t signal,Float_t signalInSigma){
	if(signalInSigma>seedSigma)
		numberOfSeeds++;
	else if(signalInSigma>hitSigma)
		numberOfHits++;
	else{
		cerr<<"No valid channel added to cluster:"<<ch<<" "<<signal<<" "<<signalInSigma<<" "<<seedSigma<<" "<<hitSigma<<(signalInSigma>hitSigma)<<endl;
		return;
	}
	cluster.push_back(make_pair(ch,signal));

}

int TCluster::size()
{
	return numberOfSeeds+numberOfHits;
}


