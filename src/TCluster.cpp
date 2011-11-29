/*
 * TCluster.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#include "../include/TCluster.hh"
ClassImp(TCluster);

TCluster::TCluster(int eventNumber,int seedSigma,int hitSigma,UInt_t nChannels) {
	// TODO Auto-generated constructor stub
	numberOfSeeds=0;
	numberOfHits=0;
	maximumSignal=0;
	charge=0;
	this->seedSigma=seedSigma;
	this->hitSigma=hitSigma;
	verbosity=0;
	revisionNumber=TCluster::TCLUSTER_REVISION();
	isChecked=false;
	isSaturated=false;
	isLumpy=false;
	isGoldenGate=false;
	hasBadChannel=false;
}

TCluster::~TCluster() {
	// TODO Auto-generated destructor stub
}

void TCluster::addChannel(int ch, Float_t signal,Float_t signalInSigma,UShort_t adcValue, bool bSaturated,bool screened){
	if(verbosity>2)cout<<"("<<ch<<"/"<<signal<<"/"<<signalInSigma<<")";
	this->isSaturated=bSaturated;
	if(signalInSigma>seedSigma)
		numberOfSeeds++;
	else if(signalInSigma>hitSigma)
		numberOfHits++;
	else{
		//cerr<<"No valid channel added to cluster:"<<ch<<" "<<signal<<" "<<signalInSigma<<" "<<seedSigma<<" "<<hitSigma<<(signalInSigma>hitSigma)<<endl;
		numberOfNoHits++;
	}
	if(signal>maximumSignal){
		maximumSignal=signal;
		maxChannel=ch;
	}
	charge+=signal;
	if(cluster.size()>0&&ch<cluster.at(0).first){
		cluster.push_back(make_pair(ch,signal));
		cluster2.push_back(make_pair(adcValue,signalInSigma));
		this->clusterChannelScreened.push_back(screened);
	}
	else{
		cluster.push_back(make_pair(ch,signal));
		cluster2.push_back(make_pair(adcValue,signalInSigma));
	}

}
Float_t TCluster::getPosition(){
	return this->getChargeWeightedMean();
	return 0;//todo;
}

UInt_t TCluster::size()
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
		if(cluster2.at(i).second>hitSigma){
			//todo . take at least second biggest hit for =charge weighted mean
			sum+=cluster.at(i).first*cluster.at(i).second;
			charged+=cluster.at(i).second;
		}
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

int TCluster::getHitSigma() const
{
    return hitSigma;
}

int TCluster::getSeedSigma() const
{
    return seedSigma;
}

void TCluster::setHitSigma(int hitSigma)
{
    this->hitSigma = hitSigma;
}

void TCluster::setSeedSigma(int seedSigma)
{
    this->seedSigma = seedSigma;
}

bool TCluster::isScreened()
{
	bool isOneChannelScreened=false;
	for(int cl;cl<clusterChannelScreened.size();cl++)
		isOneChannelScreened+=clusterChannelScreened.at(cl);
	isOneChannelScreened+=((this->getMinChannelNumber()==0)||this->getMaxChannelNumber()==nChannels-1);
}

bool TCluster::isScreened(UInt_t cl)
{
	if(cl<this->clusterChannelScreened.size())
		return clusterChannelScreened.at(cl);
	else{
		cout<<"tried to get isScreend for not valid channel from cluster:"<<cl<<endl;
		return true;
	}

}


Float_t TCluster::highest2_centroid()
{
	UInt_t maxCh = this->getMaximumChannel();
	UInt_t clPos;
	for(clPos=0;maxCh!=this->getChannel(clPos)&&clPos<size();clPos++){
	}
	if(this->getSignal(clPos-1)<getSignal(clPos+1))
		return ( getSignal(clPos)*clPos+getSignal(clPos+1)*(clPos+1))/(getSignal(clPos)+getSignal(clPos+1));
	else
		return ( getSignal(clPos)*clPos+getSignal(clPos-11)*(clPos-1))/(getSignal(clPos)+getSignal(clPos-1));

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



Float_t TCluster::getSignal(UInt_t clusterPos)
{
	if(clusterPos<size())
		return this->cluster.at(clusterPos).second;
	else return -1;
}

Float_t TCluster::getSNR(UInt_t clusterPos)
{

	if(clusterPos<size())
		return this->cluster2.at(clusterPos).second;
	else return -1;

}
Float_t TCluster::getPedestalMean(UInt_t clusterPos)
{

	if(clusterPos<size())
		return getAdcValue(clusterPos)-getSignal(clusterPos);
	else return -1;

}

Float_t TCluster::getPedestalSigma(UInt_t clusterPos)
{
	if(clusterPos<size())
		return getSignal(clusterPos)/getSNR(clusterPos);
	else return -1;

}

UShort_t TCluster::getAdcValue(UInt_t clusterPos)
{
	if(clusterPos<size())
		return this->cluster2.at(clusterPos).first;
	else return 0;
}

UInt_t TCluster::getChannel(UInt_t clusterPos)
{
	if(clusterPos<size())
		return this->cluster.at(clusterPos).first;
	else return 5000;
}
