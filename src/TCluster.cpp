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
	maximumSignal=0;
	charge=0;
	this->seedSigma=seedSigma;
	this->hitSigma=hitSigma;
	verbosity=0;
}

TCluster::~TCluster() {
	// TODO Auto-generated destructor stub
}

void TCluster::addChannel(int ch, Float_t signal,Float_t signalInSigma,bool bSaturated){
	if(verbosity>2)cout<<"("<<ch<<"/"<<signal<<"/"<<signalInSigma<<")";
	this->isSaturated=bSaturated;
	if(signalInSigma>seedSigma)
		numberOfSeeds++;
	else if(signalInSigma>hitSigma)
		numberOfHits++;
	else{
		cerr<<"No valid channel added to cluster:"<<ch<<" "<<signal<<" "<<signalInSigma<<" "<<seedSigma<<" "<<hitSigma<<(signalInSigma>hitSigma)<<endl;
		return;
	}
	if(signal>maximumSignal){
		maximumSignal=signal;
		maxChannel=ch;
	}
	charge+=signal;
	cluster.push_back(make_pair(ch,signal));

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
	return isSaturated;//todo
}
Float_t TCluster::getCharge(){
	return charge;
}
int TCluster::getMaximumChannel()
{
	return maxChannel;
}

Float_t TCluster::getChargeWeightedMean(){
	Float_t sum=0;
	Float_t charged=0;
	for(int i=0;i<this->cluster.size();i++){
		sum+=cluster.at(i).first*cluster.at(i).second;
		charged+=cluster.at(i).second;
	}

		return sum/charge;
}
void TCluster::setPositionCalulation(calculationMode_t mode){
this->mode=mode;
}

