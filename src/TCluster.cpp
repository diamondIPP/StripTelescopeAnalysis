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
	revisionNumber=TCLUSTER_REVISION;
	isChecked=false;
	isSaturated=false;
	isLumpy=false;
	isGoldenGate=false;
}

TCluster::~TCluster() {
	// TODO Auto-generated destructor stub
}

void TCluster::addChannel(int ch, Float_t signal,Float_t signalInSigma,UShort_t adcValue, bool bSaturated){
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
	if(cluster.size()>0&&ch<cluster.at(0).first){
		cluster.push_back(make_pair(ch,signal));
		cluster2.push_back(make_pair(adcValue,signalInSigma));
	}
	else{
		cluster.push_back(make_pair(ch,signal));
		cluster2.push_back(make_pair(adcValue,signalInSigma));
	}

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
bool TCluster::isLumpyCluster(){
	return isLumpy;//todo
}
bool TCluster::isGoldenGateCluster(){
	if (!isChecked)
		checkCluster();
	return this->isGoldenGate; //todo
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


void TCluster::checkCluster(){
	this->checkForGoldenGate();
	this->checkForLumpyCluster();
	isChecked=true;
}

void TCluster::checkForGoldenGate(){
	this->isGoldenGate=false;
	if(cluster.size()<=2)
		return;
	int previousSeed=-1;
	for(int i=0;i<cluster.size()&&!isGoldenGate;i++){
		if(cluster2.at(i).second>seedSigma){
			if( previousSeed!=-1 && previousSeed+1!=cluster.at(i).first )
				isGoldenGate=true;

			previousSeed=cluster.at(i).first;
		}
	}

}

void TCluster::checkForLumpyCluster(){
	this->isLumpy=false;
	if(cluster.size()<=2)
		return;//for lumpy cluster at least 3 hits are needed
	bool isfalling;
	Float_t lastSeed;
	for(int i=0;i<cluster.size()&&!isLumpy;i++){
			if(cluster2.at(i).second>seedSigma){
				if(lastSeed<cluster.at(i).second&&!isfalling){
					lastSeed=cluster.at(i).second;
				}
				if(lastSeed<cluster2.at(i).second&&isfalling)
					isLumpy=true;
				else
					isfalling=true;
			}
	}
}

bool TCluster::isSeed(UInt_t cl){
	if(cluster.size()>cl)
		return (cluster2.at(cl).second>this->seedSigma);
	return false;
}


UInt_t TCluster::getMinChannelNumber()
{
	if(size()>=0){
		return cluster.front().first;
	}
	return 0;
}

UInt_t TCluster::getMaxChannelNumber()
{

	if(size()>=0){
		return cluster.back().first;
	}
	return 0;
}
