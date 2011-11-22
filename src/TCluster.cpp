/*
 * TCluster.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#include "../include/TCluster.hh"
ClassImp(TCluster)

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
	//cluster.push_back(make_pair(ch,signal));

}
Float_t TCluster::getPosition(){
	return 0;//todo;
}

int TCluster::size()
{
	return numberOfSeeds+numberOfHits;
}

void TCluster::clear(){
	numberOfSeeds=0;
	numberOfHits=0;
	//cluster.clear();
}
bool TCluster::isLumpy(){
	return false;//todo
}
bool TCluster::isGoldenGateCluster(){
	return false; //todo
}
bool TCluster::hasSaturatedChannels(){
	return false;//todo
}
Float_t TCluster::getCharge(){
	return 0;//todo
}

void TCluster::setPositionCalulation(calculationMode_t mode){
this->mode=mode;
}

