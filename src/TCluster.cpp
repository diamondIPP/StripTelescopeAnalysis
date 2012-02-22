/*
 * TCluster.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#include "../include/TCluster.hh"
ClassImp(TCluster);

TCluster::TCluster(int nEvent,UChar_t det, int seedSigma,int hitSigma,UInt_t nChannels) {
	// TODO Auto-generated constructor stub
	clusterChannel.clear();
	clusterSignal.clear();
	clusterADC.clear();
	clusterSignalInSigma.clear();
	clusterChannelScreened.clear();
	numberOfSeeds=0;
	numberOfHits=0;
	numberOfNoHits=0;
	maximumSignal=0;
	maxChannel=0;
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
	this->det=det;
	this->eventNumber=nEvent;
	mode=highest2Centroid;
	this->nChannels=nChannels;
}

TCluster::~TCluster() {
	// TODO Auto-generated destructor stub
}


TCluster::TCluster(const TCluster& rhs){
	//TODO WIe gehe ich vor mit den ganzen TObject kram???
	verbosity=rhs.verbosity;
	if(verbosity)cout<<"copy TCluster for det  "<<(int)rhs.det<<" "<<rhs.checkClusterForSize()<<endl;
	clusterChannel.clear();
	clusterSignal.clear();
	clusterADC.clear();
	clusterSignalInSigma.clear();
	clusterChannelScreened.clear();
	numberOfSeeds=0;
	numberOfHits=0;
	numberOfNoHits=0;
	maximumSignal=0;
	maxChannel=0;
	charge=0;
	for(UInt_t i=0;i<rhs.checkClusterForSize();i++){
		UInt_t ch =rhs.clusterChannel.at(i);
		Float_t signal = rhs.clusterSignal.at(i);
		Float_t signalInSigma = rhs.clusterSignalInSigma.at(i);
		UShort_t signalAdc = rhs.clusterADC.at(i);
		clusterSignal.push_back(signal);
		clusterChannel.push_back(ch);
		clusterADC.push_back(signalAdc);
		clusterSignalInSigma.push_back(signalInSigma);
	}
	for(UInt_t i=0;i<rhs.clusterChannelScreened.size();i++)
		clusterChannelScreened.push_back(rhs.clusterChannelScreened.at(i));
	this->numberOfSeeds=rhs.numberOfSeeds;
	numberOfHits=rhs.numberOfHits;
	numberOfNoHits=rhs.numberOfNoHits;
	seedSigma=rhs.seedSigma;
	hitSigma=rhs.hitSigma;
	isSaturated=rhs.isSaturated;
	isGoldenGate=rhs.isGoldenGate;
	isLumpy=rhs.isLumpy;
	isChecked=rhs.isChecked;
	hasBadChannel=rhs.hasBadChannel;
	mode=rhs.mode;
	verbosity=rhs.verbosity;
	charge=rhs.charge;
	maximumSignal=rhs.maximumSignal;
	maxChannel=rhs.maxChannel;
	revisionNumber=rhs.revisionNumber;
	nChannels=rhs.nChannels;
	det=rhs.det;
	eventNumber=rhs.eventNumber;

}

UInt_t TCluster::checkClusterForSize() const{
	UInt_t nSignal=(clusterSignal.size());
	UInt_t nChannel=clusterChannel.size();
	UInt_t nAdc = clusterADC.size();
	UInt_t nSignalInSigma=clusterSignalInSigma.size();
	bool retVal = nSignal==nChannel;
	retVal = retVal && nChannel==nAdc;
	retVal = retVal && nAdc==nSignalInSigma;
	if(retVal==true)
		return nSignal;
	else
		return 0;
}

/**
 *
 * @param ch
 * @param signal
 * @param signalInSigma
 * @param adcValue
 * @param bSaturated
 * @param screened
 * @todo: bbSaturated not yet used....
 */
void TCluster::addChannel(UInt_t ch, Float_t signal,Float_t signalInSigma,UShort_t adcValue, bool bSaturated,bool screened){
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
	if(signalInSigma>hitSigma) charge+=signal;
	if(this->checkClusterForSize()>0&&ch<clusterChannel.front()){
		clusterChannel.push_front(ch);
		clusterSignal.push_front(signal);
		clusterSignalInSigma.push_front(signalInSigma);
		clusterADC.push_front(adcValue);
		this->clusterChannelScreened.push_front(screened);
	}
	else{
		clusterChannel.push_back(ch);
		clusterSignal.push_back(signal);
		clusterSignalInSigma.push_back(signalInSigma);
		clusterADC.push_back(adcValue);
		this->clusterChannelScreened.push_back(screened);
	}

}
Float_t TCluster::getPosition(calculationMode_t mode,TH1F *histo){
	if(mode==maxValue)
		return this->getHighestSignalChannel();
	else if(mode==chargeWeighted)
		return this->getChargeWeightedMean();
	else if(mode==highest2Centroid)
		return this->getHighest2Centroid();
	else if(mode==eta)
		return this->getEtaPostion();
	else if(mode == corEta&&histo==0)
		return this->getEtaPostion();
	else if(mode == corEta&&histo!=0)
		return this->getPositionCorEta(histo);

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
Float_t TCluster::getCharge(bool useSmallSignals){
	if(useSmallSignals)return getCharge(1000,useSmallSignals);
	return charge;
}

/**
 * @todo
 * @param nClusterEntries
 * @param useSmallSignals
 * @return
 */
Float_t TCluster::getCharge(UInt_t nClusterEntries,bool useSmallSignals){
	if(nClusterEntries==0)return 0;
	Float_t clusterCharge=getSignalOfChannel(this->getHighestSignalChannel());
	bool b2ndHighestStripIsBigger = (getHighestSignalChannel()-this->getHighest2Centroid()<0);
	for(UInt_t nCl=1;nCl<nClusterEntries&&nCl<this->checkClusterForSize();nCl++){
		if(b2ndHighestStripIsBigger){
			if(nCl%2==1){
				if (useSmallSignals||isHit(getHighestSignalChannel()+(nCl/2)+1))
					clusterCharge=clusterCharge + this->getSignalOfChannel(getHighestSignalChannel()+(nCl/2)+1);
			}
			else{
				if (useSmallSignals||isHit(getHighestSignalChannel()+(nCl/2)-1))
					clusterCharge=clusterCharge + this->getSignalOfChannel(getHighestSignalChannel()-(nCl/2));
			}
		}
		else{
			if(nCl%2==1){

				if (useSmallSignals||isHit(getHighestSignalChannel()+(nCl/2)+1))
					clusterCharge+=this->getSignalOfChannel(getHighestSignalChannel()-(nCl/2)-1);
			}
			else{
				if (useSmallSignals||isHit(getHighestSignalChannel()+(nCl/2)+1))
					clusterCharge+=this->getSignalOfChannel(getHighestSignalChannel()+(nCl/2));
			}
		}

	}
	return clusterCharge;
}

/**
 *
 * @return channel no of highest Signal
 */
UInt_t TCluster::getHighestSignalChannel()
{
	return maxChannel;
}

UInt_t TCluster::getClusterSize(){
	return this->checkClusterForSize();
}

Float_t TCluster::getChargeWeightedMean(bool useNonHits){
	Float_t sum=0;
	Float_t charged=0;
	for(UInt_t cl=0;cl<this->checkClusterForSize();cl++){
		if(useNonHits||isHit(cl)||checkClusterForSize()<4){//todo anpassen
			//todo . take at least second biggest hit for =charge weighted mean
			sum+=clusterChannel.at(cl)*clusterSignal.at(cl);//kanalnummer*signalNummer
			charged+=clusterSignal.at(cl);//signal
		}
	}

	return sum/charged;
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
	if(this->checkClusterForSize()<=2)
		return;
	int previousSeed=-1;
	for(UInt_t i=0;i<this->checkClusterForSize()&&!isGoldenGate;i++){
		if(clusterSignalInSigma.at(i)>seedSigma){
			if( previousSeed!=-1 && previousSeed+1!=(int)clusterChannel.at(i) )
				isGoldenGate=true;

			previousSeed=clusterChannel.at(i);
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
	for(UInt_t cl;cl<clusterChannelScreened.size();cl++)
		isOneChannelScreened+=clusterChannelScreened.at(cl);
	isOneChannelScreened+=((this->getSmallestChannelNumber()==0)||this->getHighestChannelNumber()==nChannels-1);
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


UInt_t TCluster::getHighestHitClusterPosition()
{
	UInt_t maxCh = this->getHighestSignalChannel();
	UInt_t clPos;
	for(clPos=0;maxCh!=this->getChannel(clPos)&&clPos<size();clPos++){
	}
	if(maxCh==getChannel(clPos))
		return clPos;
	else return 9999;
}

Float_t TCluster::getHighest2Centroid()
{
	if(getClusterSize()==0)return 0;
	UInt_t maxCh = this->getHighestSignalChannel();
	UInt_t clPos;
	for(clPos=0;maxCh!=this->getChannel(clPos)&&clPos<size();clPos++){
	}

	Float_t retVal;
	UInt_t channel=getChannel(clPos);
	Float_t signal=getSignal(clPos);
	Float_t nextChannelSignal;
	UInt_t nextChannel;
	if(this->getSignal(clPos-1)<getSignal(clPos+1)){
		nextChannel=getChannel(clPos+1);
		nextChannelSignal=getSignal(clPos+1);
//		cout<<"r"<<flush;
	}
	else{
		nextChannel=getChannel(clPos-1);
		nextChannelSignal=getSignal(clPos-1);
//		cout<<"l"<<flush;
	}
	retVal=(channel*signal+nextChannel*nextChannelSignal)/(signal+nextChannelSignal);
//	cout<<signal<<"*"<<channel<<"+"<<nextChannelSignal<<"*"<<nextChannel<<":"<<retVal<<" "<<flush;
	return retVal;
}

/**
 * todo: ueberpruefen
 */
void TCluster::checkForLumpyCluster(){
	this->isLumpy=false;
	if(this->checkClusterForSize()<=2)
		return;//for lumpy cluster at least 3 hits are needed
	bool isfalling;
	Float_t lastSeed;
	for(UInt_t i=0;i<this->checkClusterForSize()&&!isLumpy;i++){
			if(getSNR(i)>seedSigma){
				if(lastSeed<getSignal(i)&&!isfalling){
					lastSeed=getSignal(i);
				}
				if(lastSeed<getSignal(i)&&isfalling)
					isLumpy=true;
				else
					isfalling=true;
			}
	}
}

bool TCluster::isSeed(UInt_t cl){
	if(cl<checkClusterForSize())
		return (getSNR(cl)>this->seedSigma);
	return false;
}

bool TCluster::isHit(UInt_t cl){
	if(cl<checkClusterForSize())
		return (getSNR(cl)>this->hitSigma);
	return false;
}

UInt_t TCluster::getSmallestChannelNumber()
{
	if(clusterChannel.size()>0){
		return clusterChannel.front();
	}
	return 0;
}

/**
 * highest channel number of cluster...
 * goes to end of cluster and looks what is the channel no
 * @return highest channel number of cluster
 */
UInt_t TCluster::getHighestChannelNumber()
{

	if(clusterChannel.size()>0){
		return clusterChannel.back();
	}
	return 0;
}


Float_t TCluster::getHighestSignal(){
	Float_t signal=0;
	UInt_t highestSignalClusterPos;
	UInt_t highestSignalChannelPos;
	for(UInt_t clPos=0;clPos<checkClusterForSize();clPos++){
		if(getSignal(clPos)>signal){
			signal=getSignal(clPos);
			highestSignalClusterPos=clPos;
		}
	}
	highestSignalChannelPos=getChannel(highestSignalClusterPos);
	if(maximumSignal!=signal){
		cout<<"maximumSignal "<<maximumSignal<<" and highest signal "<<signal<<" does not match... Something is wrong:";
		cout<<"highestSignalChannelPos: "<<highestSignalChannelPos<<"\thighesSignalClusterPos: "<<highestSignalClusterPos<<" "<<getHighestHitClusterPosition()<<endl;

	}

	return this->maximumSignal;
}


Float_t TCluster::getSignal(UInt_t clusterPos)
{
	if(clusterPos<checkClusterForSize()){
		 Float_t signal = this->clusterSignal.at(clusterPos);
		 if(signal<0)return 0;
		 else return signal;
	}
	else {
		if(verbosity)cout<<"clusterPos "<<clusterPos<<" bigger than clusterSize"<<checkClusterForSize()<<endl;
		return 0;
	}
}

UInt_t TCluster::getClusterPosition(UInt_t channelNo){
	if(channelNo<this->getSmallestChannelNumber()&&channelNo>this->getHighestSignalChannel()){
		if(verbosity)cout<<"ChannelNo does not match..."<<endl;
		return 9999;
	}
	UInt_t clPos;
	for(clPos=0;clPos<checkClusterForSize()&&getChannel(clPos)!=channelNo;clPos++){}
	return clPos;
}

Float_t TCluster::getSignalOfChannel(UInt_t channel)
{
	if(channel<this->getSmallestChannelNumber()&&channel>this->getHighestSignalChannel()) return 0;
	UInt_t clPos;
	for(clPos=0;clPos<checkClusterForSize()&&getChannel(clPos)!=channel;clPos++){};
	if(clPos<getClusterSize()){//TODO
		Float_t signal = getSignal(clPos);
		if (signal>0)
		return signal;
	}
	return 0;
}


Float_t TCluster::getSNR(UInt_t clusterPos)
{

	if(clusterPos<checkClusterForSize())
		return this->clusterSignalInSigma.at(clusterPos);
	else return -1;

}
Float_t TCluster::getPedestalMean(UInt_t clusterPos)
{

	if(clusterPos<checkClusterForSize())
		return getAdcValue(clusterPos)-getSignal(clusterPos);
	else return -1;

}

Float_t TCluster::getPedestalSigma(UInt_t clusterPos)
{
	if(clusterPos<checkClusterForSize())
		return getSignal(clusterPos)/getSNR(clusterPos);
	else return -1;

}

UShort_t TCluster::getAdcValue(UInt_t clusterPos)
{
	if(clusterPos<checkClusterForSize())
		return this->clusterADC.at(clusterPos);
	else return 0;
}

Float_t TCluster::getEta()
{
	if (checkClusterForSize() < 2) return -1;
	UInt_t clPosHighest = getHighestHitClusterPosition();
	UInt_t clPos2ndHighest = getHighestSignalNeighbourClusterPosition(getHighestHitClusterPosition());
	UInt_t leftClPos = 0;
	UInt_t rightClPos = 0;
	if (clPosHighest < clPos2ndHighest) {
		leftClPos = clPosHighest;
		rightClPos = clPos2ndHighest;
	}
	else {
		leftClPos = clPos2ndHighest;
		rightClPos = clPosHighest;
	}
	Float_t sumSignal = (getSignal(leftClPos)+getSignal(rightClPos));
	if(sumSignal==0||getSignal(rightClPos)==0)
		return -1;
	return getSignal(rightClPos) / sumSignal;
}

Float_t TCluster::getEtaPostion(){
	if (checkClusterForSize() < 3) return -1;
	UInt_t clPosHighest = getHighestHitClusterPosition();
	UInt_t clPos2ndHighest = getHighestSignalNeighbourClusterPosition(getHighestHitClusterPosition());
	UInt_t leftClPos = 0;
	UInt_t rightClPos = 0;
	if (clPosHighest < clPos2ndHighest) {
		leftClPos = clPosHighest;
		rightClPos = clPos2ndHighest;
	}
	else {
		leftClPos = clPos2ndHighest;
		rightClPos = clPosHighest;
	}
	Float_t eta = getEta();
	if(eta==0)return -9999;
	return eta+leftClPos;
}

Float_t TCluster::getPositionCorEta(TH1F* histo){
	if (checkClusterForSize() < 3) return -1;
	if(histo==0) return -1;
	UInt_t clPosHighest = getHighestHitClusterPosition();
	UInt_t clPos2ndHighest = getHighestSignalNeighbourClusterPosition(getHighestHitClusterPosition());
	UInt_t leftClPos = 0;
	UInt_t rightClPos = 0;

	Float_t eta = getEta();
	if(verbosity)cout<<"get Position Cor Eta: "<<eta<<" ";
	if (clPosHighest < clPos2ndHighest) {
		leftClPos = clPosHighest;
		rightClPos = clPos2ndHighest;
		if(verbosity)cout<<"leftHighest "<<flush;
	}
	else {
		leftClPos = clPos2ndHighest;
		rightClPos = clPosHighest;
		if(verbosity)cout<<"rightHighest "<<flush;
	}
	if(eta<=0)
		return -1;
	Float_t corEta= getValueOfHisto(eta,histo);
	UInt_t leftChannelNo= getChannel(leftClPos);
	if(verbosity)	cout<<leftChannelNo<<" + "<<corEta<<" = "<<leftChannelNo+corEta;
	return leftChannelNo+corEta;
}
/**
 * todo: return value 5000 check if that makes sense
 * @param clusterPos
 * @return
 */
UInt_t TCluster::getChannel(UInt_t clusterPos)
{
	if(clusterPos<checkClusterForSize())
		return this->clusterChannel.at(clusterPos);
	else return 5000;
}


UInt_t TCluster::getHighestSignalNeighbourChannel(UInt_t channelNo)
{
	return getChannel(getHighestSignalNeighbourClusterPosition(getClusterPosition(channelNo)));
}

UInt_t TCluster::getHighestSignalNeighbourClusterPosition(UInt_t clPos)
{
	if (clPos>=checkClusterForSize() || clPos<0 || checkClusterForSize()<2) return 9999;
	if (getSignal(clPos-1) < getSignal(clPos+1)) return clPos+1;
	if (getSignal(clPos-1) > 0) return clPos-1;
	return 9999;
}

/**
 * small function to Intend cout-output;
 * @input level level of intention
 * @return string with level tabs
 * @todo move to "string stuff" class
 */
string TCluster::Intent(UInt_t level){
	stringstream output;
	output.str("");
	for(UInt_t i=0;i<level;i++){
		output<<"  ";
	}
	return output.str();
}
Float_t TCluster::getValueOfHisto(Float_t x, TH1F* histo){
	if(histo->IsZombie()){
		cerr<<"Histo is Zombie!"<<endl;
		return -999;
	}
	Float_t xmin = histo->GetXaxis()->GetXmin();
	Float_t xmax = histo->GetXaxis()->GetXmax();
	if(xmin<=x&&x<=xmax){
		Int_t bin = histo->FindBin(x);
		return histo->GetBinContent(bin);
	}
	cerr<<" x = "<<x<<" not in range ["<<xmin<<","<<xmax<<"]"<<endl;
	return -999;

}

void TCluster::Print(UInt_t level){
	cout<<Intent(level)<<"Cluster of Event "<<flush;
	cout<<eventNumber<<" in detector"<<(int)det<<" with "<<size()<<" Cluster entries"<<flush;
	for(UInt_t cl=0;cl<checkClusterForSize();cl++){
		if(this->isSeed(cl))
			cout<<"\t{"<<this->getChannel(cl)<<"|"<<this->getAdcValue(cl)<<"|"<<this->getSignal(cl)<<"|"<<this->getSNR(cl)<<"}"<<flush;
		else if(this->isHit(cl))
					cout<<"\t("<<this->getChannel(cl)<<"|"<<this->getAdcValue(cl)<<"|"<<this->getSignal(cl)<<"|"<<this->getSNR(cl)<<")"<<flush;
		else
					cout<<"\t["<<this->getChannel(cl)<<"|"<<this->getAdcValue(cl)<<"|"<<this->getSignal(cl)<<"|"<<this->getSNR(cl)<<"]"<<flush;
	}
	cout<<"\t||"<<this->getHighestSignalChannel()<<" "<<flush<<this->getHighest2Centroid()<<" "<<this->getChargeWeightedMean();
	cout<<endl;
}
