/*
 * TCluster.cpp
 *
 *  Created on: 21.11.2011
 *      Author: bachmair
 */

#include "../include/TCluster.hh"
ClassImp(TCluster);

TCluster::TCluster(int nEvent,UChar_t det, int seedSigma,int hitSigma,UInt_t nChannels,Float_t cmNoise) {
	initialiseNewCluster();
	this->seedSigma=seedSigma;
	this->hitSigma=hitSigma;
	this->det=det;
	this->eventNumber=nEvent;
	this->nChannels=nChannels;
	this->cmNoise = cmNoise;
	if (verbosity>2) cout<<"new TCluster of event "<<nEvent<<" of detector "<<(int)det<<" Common mode Noise "<<cmNoise<<endl;
}

void TCluster::initialiseNewCluster(){
    verbosity = 0;
	clusterChannel.clear();
	clusterSignal.clear();
	clusterADC.clear();
	clusterSignalInSigma.clear();
	clusterChannelScreened.clear();
	clusterPedMean.clear();
	clusterPedSigma.clear();
	clusterPedMeanCMN.clear();
	clusterPedSigmaCMN.clear();
    numberOfSeeds = 0;
    numberOfHits = 0;
    seedSigma = 10;
    hitSigma = 7;
    isSaturated = false;
    isGoldenGate = false;
    isLumpy = false;
    maximumSignal = 0;
    charge = 0;
    revisionNumber=TCLUSTER_REVISION();
    isChecked = false;
    hasBadChannel=false;
    numberOfNoHits=0;
    nChannels=256;
    isTransparentCluster = -1;
    cmNoise = 0;
    mode=highest2Centroid;
    transparentClusterSize =0;
}

TCluster::TCluster(const TCluster& rhs){
	initialiseNewCluster();
	//TODO WIe gehe ich vor mit den ganzen TObject kram???
	verbosity=rhs.verbosity;
	if(verbosity)cout<<"copy TCluster for det  "<<(int)rhs.det<<" "<<rhs.checkClusterForSize()<<endl;
	cmNoise=rhs.cmNoise;
	for(UInt_t i=0;i<rhs.checkClusterForSize();i++){
		clusterSignal.push_back(rhs.clusterSignal.at(i));
		clusterChannel.push_back(rhs.clusterChannel.at(i));
		clusterADC.push_back( rhs.clusterADC.at(i));
		clusterSignalInSigma.push_back( rhs.clusterSignalInSigma.at(i));
		clusterPedMean.push_back(rhs.clusterPedMean.at(i));
		clusterPedSigma.push_back(rhs.clusterPedSigma.at(i));
		clusterPedMeanCMN.push_back(rhs.clusterPedMeanCMN.at(i));
		clusterPedSigmaCMN.push_back(rhs.clusterPedSigmaCMN.at(i));
	}
	for(UInt_t i=0;i<rhs.clusterChannelScreened.size();i++)
		clusterChannelScreened.push_back(rhs.clusterChannelScreened[i]);
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
	isTransparentCluster=rhs.isTransparentCluster;
	transparentClusterSize=rhs.transparentClusterSize;
	//cout<<"transparentClusterSize"<<transparentClusterSize<<"/"<<rhs.transparentClusterSize<<endl;
}

/**
 * the class assignment Function
 * @param src
 * @return
 */
TCluster & TCluster::operator =(const TCluster & src)
{
	//To clear and empty all vectors.
	clusterChannel.clear();clusterChannel.empty();
	clusterSignal.clear();clusterSignal.empty();
	clusterADC.clear();clusterADC.empty();
	clusterSignalInSigma.clear();clusterSignalInSigma.empty();
	clusterPedMean.clear();clusterPedMean.empty();
	clusterPedSigma.clear();clusterPedSigma.empty();
	clusterPedMeanCMN.clear();clusterPedMeanCMN.empty();
	clusterPedSigmaCMN.clear();clusterPedSigmaCMN.empty();
	clusterChannelScreened.clear();clusterChannelScreened.empty();

	hitSigma = src.hitSigma;
	seedSigma= src.seedSigma;
	isSaturated= src.isSaturated;
	isLumpy=src.isLumpy;
	isGoldenGate=src.isGoldenGate;
	isChecked=src.isChecked;
	hasBadChannel=src.hasBadChannel;
	mode=src.mode;
	numberOfHits=src.numberOfHits;
	numberOfNoHits=src.numberOfNoHits;
	numberOfSeeds=src.numberOfSeeds;
	verbosity=src.verbosity;
	charge=src.charge;
	//printf("srcCharge: %f, CopiedCharge: %f \n",src.charge,charge);
	maxChannel=src.maxChannel;
	maximumSignal=src.maximumSignal;
	revisionNumber=src.revisionNumber;
	nChannels=src.nChannels;
	det=src.det;
	eventNumber=src.eventNumber;
	//printf("srceventNumber: %f, CopiedeventNumber: %f \n",src.eventNumber,eventNumber);
	cmNoise= src.cmNoise;
	isTransparentCluster = src.isTransparentCluster;
	//printf("srcisTransparentCluster: %f, CopiedisTransparentCluster: %f \n",src.isTransparentCluster,isTransparentCluster);
	transparentClusterSize = src.transparentClusterSize;
	//cout<<"transparentClusterSize: "<<transparentClusterSize<<"/"<<src.transparentClusterSize<<endl;
	//printf("srcClusterSize: %i, CopiedClusterSize: %i \n",src.transparentClusterSize,transparentClusterSize);
	for(UInt_t i=0;i<src.clusterChannel.size();i++)
		clusterChannel.push_back(src.clusterChannel.at(i));
	for(UInt_t i=0;i<src.clusterSignal.size();i++)
		clusterSignal.push_back(src.clusterSignal.at(i));
	for(UInt_t i=0;i<src.clusterADC.size();i++)
		clusterADC.push_back(src.clusterADC.at(i));
	for(UInt_t i=0;i<src.clusterSignalInSigma.size();i++)
		clusterSignalInSigma.push_back(src.clusterSignalInSigma.at(i));
	for(UInt_t i=0;i<src.clusterPedMean.size();i++)
		clusterPedMean.push_back(src.clusterPedMean.at(i));
	for(UInt_t i=0; i<src.clusterPedSigma.size();i++)
		clusterPedSigma.push_back(src.clusterPedSigma.at(i));
	for(UInt_t i=0;i<src.clusterPedMeanCMN.size();i++)
		clusterPedMeanCMN.push_back(src.clusterPedMeanCMN.at(i));
	for(UInt_t i=0; i<src.clusterPedSigmaCMN.size();i++)
		clusterPedSigmaCMN.push_back(src.clusterPedSigmaCMN.at(i));
	for(UInt_t i=0;i<src.clusterChannelScreened.size();i++)
		clusterChannelScreened.push_back(src.clusterChannelScreened.at(i));
	return *this;
}


TCluster::~TCluster() {
	while (clusterSignal.size())	clusterSignal.pop_back();
	while (clusterChannel.size()) clusterChannel.pop_back();
	while (clusterADC.size()) clusterADC.pop_back();
	while (clusterSignalInSigma.size()) clusterSignalInSigma.pop_back();
	while (clusterPedMean.size()) clusterPedMean.pop_back();
	while (clusterPedSigma.size()) clusterPedSigma.pop_back();
	while (clusterPedMeanCMN.size()) clusterPedMeanCMN.pop_back();
	while (clusterPedSigmaCMN.size()) clusterPedSigmaCMN.pop_back();
}

void TCluster::SetTransparentCluster(Float_t startChannel) {
	if(startChannel>=0&&startChannel< nChannels)
		isTransparentCluster=startChannel;
	if(startChannel==-1)
		isTransparentCluster=startChannel;
}
void TCluster::SetTransparentClusterSize(UInt_t size){
    if(size>0) transparentClusterSize=TMath::Min(size,checkClusterForSize());
//    cout<<"TCluster:SetTransparentClusterSize("<<size<<"): - "<<transparentClusterSize<<endl;
//    cout<<"TCluster:SetTransparentClusterSize("<<size<<"): - "<<transparentClusterSize<<endl;
//    char t;
//    cin>>t;
};

Int_t TCluster::getTransparentClusterPosition(UInt_t clusterNo) {
	Int_t startChannel = Int_t(isTransparentCluster+.5);
	Int_t clStart = this->getClusterPosition(startChannel);
	Int_t dir;
	if (startChannel-isTransparentCluster<0)
	    dir = -1;
	else
	    dir = 1;
	if (clusterNo%2==0) dir *= -1;
	Int_t dif = (Int_t(clusterNo)+1)/2;
//	if(isTransparentCluster-startChannel>0){
//		if(clusterNo%2==0)
//			dir = -1;
//	}
//	else{
//		if(clusterNo%2==1)
//			dir = -1;
//	}
	Int_t clusPos = clStart + dir *dif;
    //cout<<"[TCluster::getTransparentClusterPosition] "<<startChannel<<" "<<clStart<<" "<<dir<<" "<<dif<<" "<<clusPos<<endl;
	return clusPos;
}

/**
 * checks if all deques with the data of the cluster do have the same size
 * @return the size of the cluster if all sizes are the same, else 0
 */
UInt_t TCluster::checkClusterForSize() const{
	UInt_t nChannel=clusterChannel.size();
	UInt_t nAdc = clusterADC.size();
	UInt_t nPedMean = clusterPedMean.size();
	UInt_t nPedMeanCMN = clusterPedMeanCMN.size();
	UInt_t nPedSigma = clusterPedSigma.size();
	UInt_t nPedSigmaCMN = clusterPedSigmaCMN.size();

	bool retVal = (nChannel==nAdc);
	retVal = retVal && nAdc==nPedMean;
	retVal = retVal && nPedMean==nPedMeanCMN;
	retVal = retVal && nPedMeanCMN==nPedSigma;
	retVal = retVal && nPedSigma==nPedSigmaCMN;
	if(retVal==true)
		return nAdc;
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

void TCluster::addChannel(UInt_t ch, Float_t pedMean, Float_t pedSigma, Float_t pedMeanCMN, Float_t pedSigmaCMN, Int_t adcValue, bool bSaturated,bool isScreened){
	//void TCluster::addChannel(UInt_t ch, Float_t signal,Float_t signalInSigma,UShort_t adcValue, bool bSaturated,bool screened){
	if(ch>TPlaneProperties::getNChannels((UInt_t)this->det)||ch!=ch){
		cerr<<"The channel which should be added is not VALID! you cannot add a channel "<<ch<<" at det "<<det<<endl;
		return;
	}
	if (pedMean!=pedMean||pedSigma!=pedSigma||pedMeanCMN!=pedMeanCMN||pedSigmaCMN!=pedSigmaCMN||adcValue!=adcValue){
		cerr<<"Cannot add channel: one of the variables is invalid"<<endl;
		return;
	}
	Float_t signal = adcValue - pedMean;
	Float_t signalInSigma= signal/pedSigma;
	Float_t signalCMN = adcValue - pedMeanCMN - this->cmNoise;
	Float_t signalInSigmaCMN= signal/pedSigma;
	if(verbosity>2)cout<<"("<<ch<<"/"<<signal<<"/"<<signalInSigma<<")";
	this->isSaturated=this->isSaturated||bSaturated;
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
		clusterSignalCMN.push_back(signalCMN);
		clusterSignalInSigmaCMN.push_back(signalInSigmaCMN);
		clusterADC.push_front(adcValue);
		this->clusterChannelScreened.push_front(isScreened);
		clusterPedMean.push_front(pedMean);
		clusterPedMeanCMN.push_front(pedMeanCMN);
		clusterPedSigma.push_front(pedSigma);
		clusterPedSigmaCMN.push_front(pedSigmaCMN);
	}
	else{
		clusterChannel.push_back(ch);
		clusterSignal.push_back(signal);
		clusterSignalInSigma.push_back(signalInSigma);
		clusterSignalCMN.push_back(signalCMN);
		clusterSignalInSigmaCMN.push_back(signalInSigmaCMN);
		clusterADC.push_back(adcValue);
		clusterPedMean.push_back(pedMean);
		clusterPedMeanCMN.push_back(pedMeanCMN);
		clusterPedSigma.push_back(pedSigma);
		clusterPedSigmaCMN.push_back(pedSigmaCMN);
		this->clusterChannelScreened.push_back(isScreened);
	}
	if(verbosity>2)cout<<flush;

}

/**
 * Position of Cluster in channels calculated by calculation mode
 * @param mode
 * @param histo needed for eta corrected calculation
 * @return
 */
Float_t TCluster::getPosition(bool cmnCorrected, calculationMode_t mode,TH1F *histo){
	if(mode == TCluster::corEta && histo==0)
		cout<<"TCluster::getPosition::histo == 0"<<endl;
	if(mode==maxValue)
		return this->getHighestSignalChannel();
	else if(mode==chargeWeighted)
		return this->getChargeWeightedMean(cmnCorrected);
	else if (mode == highest2CentroidNoSmallHits)
		return this->getHighest2Centroid(cmnCorrected,false);
	else if(mode==highest2Centroid)
		return this->getHighest2Centroid(cmnCorrected,true);
	else if(mode==eta)
		return this->getEtaPostion(cmnCorrected);
	else if(mode == corEta&&histo==0){
		if(verbosity)cerr<<"mode = cor Eta, but histo =0, "<<endl;
		//		this->Print();
		return this->getEtaPostion(cmnCorrected);
	}
	else if(mode == corEta&&histo!=0)
		return this->getPositionCorEta(cmnCorrected,histo);

	return 0;//todo;
}

UInt_t TCluster::size()
{
    if (IsTransparentCluster())
        return GetMaxTransparentClusterSize();
	return numberOfSeeds+numberOfHits;
}

UInt_t TCluster::seedSize(){
	return numberOfSeeds;
}
void TCluster::clear(){
	numberOfSeeds=0;
	numberOfHits=0;
	clusterChannel.clear();
	clusterSignal.clear();
	clusterADC.clear();
	clusterSignalInSigma.clear();
	clusterChannelScreened.clear();
	clusterPedMean.clear();
	clusterPedSigma.clear();
	clusterPedMeanCMN.clear();
	clusterPedSigmaCMN.clear();
	//cluster.clear();
}
bool TCluster::isLumpyCluster(){
	if (!isChecked)
		checkCluster();
	return isLumpy;
}

bool TCluster::isGoldenGateCluster(){
	if (!isChecked)
		checkCluster();
	return this->isGoldenGate;
}

bool TCluster::hasSaturatedChannels(){
	if(IsTransparentCluster()){
		for(UInt_t i=0; i<GetTransparentClusterSize();i++)
				if (getAdcValue(getTransparentClusterPosition(i))>=TPlaneProperties::getMaxSignalHeight(this->det)) return true;
	}
	else{
		for(UInt_t cl=0;cl<this->clusterADC.size();cl++){
			if (getAdcValue(cl)>=TPlaneProperties::getMaxSignalHeight(this->det)) return true;
		}
	}
	return false;

	return isSaturated;//todo
}

TCluster TCluster::getCrossTalkCorrectedCluster(Float_t alpha){
	TCluster newClus = TCluster(this->eventNumber,this->det,this->seedSigma,this->hitSigma,this->nChannels,this->cmNoise);
	newClus.clear();//(clus.getEventNumber(),clus.getD);
	UInt_t det = this->getDetector();
	UInt_t clSize = this->getClusterSize();
	Float_t S_i = 0;
	bool do_cmn = false;
	if (det == TPlaneProperties::isDiamondDetector(det))
	    do_cmn = true;
	TString str = "";
	Int_t direction = copysign(1.0,alpha);
	alpha = copysign(alpha,1.0);
	str += TString::Format("\ngetCrossTalkCorrectedCluster for det %d: %5.2f %%\n",det,alpha*100);
    bool bChanged = false;
	for(UInt_t cl = 0; cl < clSize;cl++){
		UInt_t clPos = cl;
		if(direction == -1)
			clPos = clSize - cl-1;
		Int_t adc = this->getAdcValue(clPos);
		Float_t ped = this->getPedestalMean(clPos,do_cmn);
		Float_t measured_signal = adc-ped;

		UInt_t channel = this->getChannel(clPos);
		str+= TString::Format("\t*%2d - %3d --> %4d - %6.1f = %6.1f, %6.1f", cl,channel,(Int_t)adc,ped,measured_signal,S_i);
        S_i = (measured_signal - S_i * alpha)/(1-alpha);
        Int_t old_adc = adc;
        adc = (Int_t)(S_i+ped+0.5);
        str+=  TString::Format(" ==> %6.1f / %4d  - %1d\n", S_i,adc,(old_adc-adc));
		bool isSaturated = this->getAdcValue(clPos)>=TPlaneProperties::getMaxSignalHeight(det);
		newClus.addChannel(channel,this->getPedestalMean(clPos),this->getPedestalSigma(clPos),
		        this->getPedestalMean(clPos,true),this->getPedestalSigma(clPos,true),adc,
				isSaturated,this->isScreened(clPos));
        if (old_adc!=adc) bChanged = true;
	}
    if (bChanged&&false){
	    if (alpha) cout<<str<<endl;
	    if (alpha) cout<<"OLD: "<<endl;
	    if (alpha) this->Print(1);
	    if (alpha) cout<<"NEW: "<<endl;
	    if (alpha) newClus.Print(1);
    }
	return newClus;
}

bool TCluster::IsValidTransparentClusterPosition(UInt_t clusPos){
    Int_t cp1 = GetTransparentClusterSize()-1;
    Int_t cp2 = GetTransparentClusterSize()-2;

	Int_t cl1 = this->getTransparentClusterPosition(cp1);
	Int_t cl2 = this->getTransparentClusterPosition(cp2);
	if(cl1<0)
		cl1 = 0;
	if(cl2<0)
		cl2 = 0;
	if (cl2 > this->clusterADC.size())
		cl2 = clusterADC.size();
	if (cl1 > this->clusterADC.size())
			cl1 = clusterADC.size();
	if(cl2<cl1){
		Int_t cl3 = cl2;
		cl2 = cl1;
		cl1 = cl3;
	}
	if(cl1<=clusPos&&clusPos<=cl2)
		return true;
    if (false && clusPos< clusterADC.size()){
        cout<<"TCluster::IsValidTransparentClusterPosition("<<clusPos<<"):"<<GetTransparentClusterSize()<<"/"<<flush;
        cout<<"cl: ["<<cl1<<"-"<<clusPos<<"-"<<cl2<<"]"<<flush;
        cout<<"                     "<<this->getTransparentClusterPosition(cp1)<<"/"<<
            this->getTransparentClusterPosition(cp2) <<"   "<<cp1<<"//"<<cp2<<endl;
        for (UInt_t i = 0; i < clusterADC.size();i++)
            cout<<" "<<setw(2)<<i<<": "<<this->getTransparentClusterPosition(i)<<" "<<
                getAdcValue(this->getTransparentClusterPosition(i))<<endl;
    }
    //this->Print();
	return false;
}

TCluster::direction_t TCluster::getDirection(UInt_t clusterPosition){

	if(IsTransparentCluster()){
		if(!IsValidTransparentClusterPosition(clusterPosition-1))
			return right;
		if(!IsValidTransparentClusterPosition(clusterPosition+1))
			return left;
		if(getSignal(clusterPosition-1)<getSignal(clusterPosition+1))
			return right;
		else
			return left;
	}
	else{
		if(getSignal(clusterPosition-1)<getSignal(clusterPosition+1))
			return right;
		else
			return left;
	}
}


Float_t TCluster::getTransparentCharge(UInt_t clusters, bool cmnCorrected,
		bool useSmallSignals) {
	if(!IsTransparentCluster())
		return getCharge(clusters,cmnCorrected,useSmallSignals);
//	Int_t channel = (Int_t)(isTransparentCluster+.5);
//	direction_t direction = isTransparentCluster-channel<0?left:right;
	UInt_t startingClusterPos = getHighestHitClusterPosition();
	direction_t direction = getDirection(startingClusterPos);
	if(clusters> GetTransparentClusterSize())
		clusters = GetTransparentClusterSize();
	if(startingClusterPos < checkClusterForSize())
		return getChargeStartingAt(clusters,startingClusterPos,direction,false,cmnCorrected,useSmallSignals);
	return 0;

}

Float_t TCluster::getChargeStartingAt(UInt_t nChannels, UInt_t startingClusterPos,direction_t direction,
        bool noNegativeCharge, bool useCMcorrection, bool useSmallSignals) {
//	if(verbosity>3)
	UInt_t size = checkClusterForSize();
	nChannels = TMath::Min(nChannels,2*size);
	if(IsTransparentCluster()){
		nChannels = TMath::Min(nChannels,GetTransparentClusterSize());
	}
	if(verbosity>3)cout<<" "<<nChannels<<endl;
	if(startingClusterPos>=size)
		return 0;
	Float_t charge;
	Float_t clusterCharge = 0;
	Int_t channel;
	UInt_t nUsedChannels=0;
	bool useEvenClusterPos =  true;
	bool useOddClusterPos = true;
//	getSignal(startingClusterPos,useCMcorrection);

	for(UInt_t nCl = 0; nUsedChannels < nChannels&& nCl<2*size; nCl++){
		UInt_t clPos = startingClusterPos;
		Int_t distance = ((nCl+1)/2);
		Int_t sign = 1;
		if((direction==left&&nCl%2==1)||(direction == right && nCl%2 == 0))
			sign =-1;
		clPos += sign * distance;
		channel = getChannel(clPos);
//		if(nChannels ==7 && eventNumber == 40858)cout<<"\n "<<nUsedChannels<<"/"<<nChannels<<"\t"<<nCl<<" "<<clPos<<":\t"<<channel<<"\t";
		if(IsTransparentCluster()&&!IsValidTransparentClusterPosition(clPos))
			continue;
		if(clPos>= size)
			continue;
		charge = 0;
		if (useSmallSignals|| isHit(clPos))
			charge = getSignal(clPos,useCMcorrection);
		if(charge < 0 && noNegativeCharge){
		    if(clPos%2==0)
		        useEvenClusterPos = false;
		    else
		        useOddClusterPos = false;
		    charge = 0;
		}
		if(noNegativeCharge && !useOddClusterPos && clPos%2==1)
		    charge = 0;
		if(noNegativeCharge && !useEvenClusterPos && clPos%2==0)
		    charge =0;
		clusterCharge+=charge;
//		cout<<" "<<charge<<" "<<clusterCharge;
		nUsedChannels ++;
//		if(verbosity>3)cout<<TString::Format("%2d %2d %2d %2d",nCl,sign,distance,clPos)<<"\t"<<channel<<"-->"<<charge<<" "<<clusterCharge<<endl;
	}
//	cout<<endl;
	return clusterCharge;
}

Float_t TCluster::getCharge(bool cmnCorrected,bool useSmallSignals){
//	cout<<"getCHarge"<<endl;
	//	if(useSmallSignals)
	return getCharge(1000,cmnCorrected,useSmallSignals);
	//	return charge;
}

Float_t TCluster::getPositiveCharge(bool cmnCorrected,bool useSmallSignals){
//  cout<<"getCHarge"<<endl;
    //  if(useSmallSignals)
    return getPositiveCharge(1000,cmnCorrected,useSmallSignals);
    //  return charge;
}
/**
 * @param nClusterEntries
 * @param useSmallSignals
 * @return
 */
Float_t TCluster::getCharge(UInt_t nClusterEntries,bool cmnCorrected,bool useSmallSignals){
	if(nClusterEntries==0)return 0;
	Int_t clusPosStart ;
	direction_t dir;
	if (!IsTransparentCluster()){
	    clusPosStart= getClusterPosition(getHighestSignalChannel());
	    dir = getSignal(clusPosStart-1,cmnCorrected)>getSignal(clusPosStart+1,cmnCorrected)?left:right;
	}
	else{
	    clusPosStart= getClusterPosition((int)(isTransparentCluster+.5));
	    dir = isTransparentCluster-(int)(isTransparentCluster+.5)>0?right:left;
	}
	Float_t clusterCharge = getChargeStartingAt(nClusterEntries,clusPosStart,dir,false,cmnCorrected,useSmallSignals);
	return clusterCharge;
}


//
UInt_t TCluster::GetHighestSignalChannelTransparentCluster(){
////	if(!IsTransparentCluster())
////		return getHighestChannelNumber();
//	if(checkClusterForSize()<10)
//		return -999;
//	cout<<checkClusterForSize()<<endl;
//	Print(1);
//	cout<<endl;
	UInt_t ch = isTransparentCluster+.5;
	UInt_t clPosStart = getClusterPosition(ch);
	Int_t maxChannel = -1;
	Float_t maxSignal = -1e9;
	Int_t dir = ch-isTransparentCluster>0?-1:+1;
	UInt_t len = TMath::Min(GetTransparentClusterSize(),2*this->checkClusterForSize());
//	cout<<"GetHighestSignalChannelTransparentCluster: "<< ch<<" "<<clPosStart<< " "<< dir<<" "<<len<<endl;
	for (UInt_t i = 0; i< len; i++){
		UInt_t clPos = getTransparentClusterPosition(i);
//		Int_t j = (i+1)/2;
//		dir *=-1;
//		UInt_t clPos = clPosStart+dir*j;
//		cout<< i <<" "<<j<<" "<<dir<<" "<<clPos<<endl;
		if(getSignal(clPos)>maxSignal){
			maxSignal = getSignal(clPos);
			maxChannel = getChannel(clPos);
		}

	}
	return maxChannel;
	return -9999;
}
/**
 *
 * @return channel no of highest Signal
 */
UInt_t TCluster::getHighestSignalChannel()
{
	if(IsTransparentCluster()){
		return GetHighestSignalChannelTransparentCluster();
	}
	return maxChannel;
}


void TCluster::UpdateHighestSignalChannel() {

	Float_t maxSignal = getSignal(0);
	Int_t maxChannel = getChannel(0);
	for(UInt_t cl=1;cl<this->checkClusterForSize();cl++){
		if(getSignal(cl)>maxSignal){
			maxChannel = getChannel(cl);
			maxSignal=getSignal(cl);
		}
	}
	this->maxChannel=maxChannel;
}

UInt_t TCluster::getClusterSize(){
	if(IsTransparentCluster())
		return TMath::Min(GetTransparentClusterSize(),checkClusterForSize());
	return this->checkClusterForSize();
}

Float_t TCluster::getChargeWeightedMean(bool cmnCorrected, bool useNonHits){
	Float_t sum=0;
	Float_t charged=0;
	for(UInt_t cl=0;cl<this->checkClusterForSize();cl++){
		if(IsTransparentCluster()&&!IsValidTransparentClusterPosition(cl))
			continue;
		if(useNonHits||isHit(cl)||checkClusterForSize()<4){//todo anpassen
			//todo . take at least second biggest hit for =charge weighted mean
			sum+=getChannel(cl)*getSignal(cl,cmnCorrected);//kanalnummer*signalNummer
			charged+=getSignal(cl,cmnCorrected);;//signal
		}
	}
	if (charged>0)
		return sum/charged;
	return -1;
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
	if(verbosity>10)cout<< "check for golden gate cluster!: "<<endl;
	this->isGoldenGate=false;
	if(this->checkClusterForSize()<=2)
		return;
	if(verbosity>10)cout<<"clSize:"<<checkClusterForSize()<<" "<<clusterChannel.size()<<flush;
	int previousSeed=-1;
	for(UInt_t i=0;i<this->checkClusterForSize()&&!isGoldenGate;i++){
		if(getSNR(i)>seedSigma){
			if( previousSeed!=-1 && previousSeed+1!=(int)getChannel(i) )
				isGoldenGate=true;

			previousSeed=(int)getChannel(i);
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
	for(UInt_t cl = 0;cl<clusterChannelScreened.size();cl++)
		isOneChannelScreened+=isScreened(cl);
	isOneChannelScreened+=((this->getSmallestChannelNumber()==0)||this->getHighestChannelNumber()==nChannels-1);
	return isOneChannelScreened;
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
UInt_t TCluster::getHighestHitClusterPositionTransparentCluster(){
		UInt_t ch = isTransparentCluster+.5;
		Int_t maxClPos = -1;
		Float_t maxSignal = -1e9;
		UInt_t len = TMath::Min(GetTransparentClusterSize(),2*this->checkClusterForSize());
		for (UInt_t i = 0; i< len; i++){
			UInt_t clPos = getTransparentClusterPosition(i);
			if(getSignal(clPos)>maxSignal){
				maxSignal = getSignal(clPos);
				maxClPos = clPos;
			}

		}
		return maxClPos;
		return -9999;
}

UInt_t TCluster::getHighestHitClusterPosition()
{
	if(IsTransparentCluster()){
		return getHighestHitClusterPositionTransparentCluster();
	}
	UInt_t maxCh = this->getHighestSignalChannel();

	UInt_t clPos;
	for(clPos=0;maxCh!=this->getChannel(clPos)&&clPos<size();clPos++){
	}
	if(maxCh==getChannel(clPos))
		return clPos;
	else return 9999;
}
/** should be the same as getEtaPosition
 *
 * @return
 */
Float_t TCluster::getHighest2Centroid(bool cmnCorrected, bool useSmallSignals)
{
	if(getClusterSize()==0)return 0;
	UInt_t maxCh = this->getHighestSignalChannel();
	UInt_t clPos = getClusterPosition(maxCh);

	Float_t retVal;
	UInt_t channel=getChannel(clPos);
	Float_t signal=getSignal(clPos,cmnCorrected);
	Float_t nextChannelSignal;
	UInt_t nextChannel;
	Float_t signalLeft = getSignal(clPos-1,cmnCorrected);
	Float_t signalRight = getSignal(clPos+1,cmnCorrected);
	if(IsTransparentCluster()){
		if(!IsValidTransparentClusterPosition(clPos-1))
			signalLeft=-1;
		if(!IsValidTransparentClusterPosition(clPos+1))
			signalRight=-1;
	}
	if (signalLeft<0&&signalRight<0)
			return maxCh;
	if(!useSmallSignals && !isHit(clPos-1) && !isHit(clPos+1))
		return maxCh;
	if(signalLeft<signalRight){
		nextChannel=getChannel(clPos+1);
		nextChannelSignal=signalRight;
		//		cout<<"r"<<flush;
	}
	else{
		nextChannel=getChannel(clPos-1);
		nextChannelSignal=signalLeft;
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
	if(verbosity>10)cout<<"check is lumpy:"<<endl;
	this->isLumpy=false;
	UInt_t clSize = checkClusterForSize();
	if(this->checkClusterForSize()<=2)
		return;//for lumpy cluster at least 3 hits are needed
	bool isfalling;
	Float_t lastSeed = -99999.;
	if(verbosity>10)cout<<" Do Loop: "<<endl;
	for(UInt_t i=0;i<clSize&&!isLumpy;i++){
		float signal = getSignal(i);
		float snr = getSNR(i);
		if(verbosity>10)cout<<" "<<i<<flush;
		if(snr>seedSigma){
			if(verbosity>10)cout<<" - found seed with SNR "<<snr<<flush;
			if(lastSeed<signal&&!isfalling){
				lastSeed=getSignal(i);
				if(verbosity>10)cout<<", found new las seed with signal "<< getSignal(i)<<flush;
			}
			if(lastSeed<signal&&isfalling){
				isLumpy=true;
				if(verbosity>10)cout<<", is LUMPY"<<flush;
			}
			else{
				isfalling=true;
				if(verbosity>10)cout<<", is falling"<<flush;
			}
		}
		if(verbosity>10)cout<<endl;
	}
	if(verbosity>10)cout<<"isLumpyCluster: "<<isLumpy<<endl;
}

/**
 * checks if an entry in the cluster is above the seed SNR
 * @param cl
 * @return
 */
bool TCluster::isSeed(UInt_t cl){
	if(cl<checkClusterForSize())
		return (getSNR(cl)>this->seedSigma);
	return false;
}
/**
 * checks if an entry in the cluster is above the hit SNR
 * Could be either a hit or a seed
 * @param cl
 * @return
 */
bool TCluster::isHit(UInt_t cl){
	if(cl<checkClusterForSize())
		return (getSNR(cl)>this->hitSigma);
	return false;
}

UInt_t TCluster::getSmallestChannelNumber()
{
	if(clusterChannel.size()==0)
		return 0;
//	if(IsTransparentCluster()){
//		Int_t channel1 = getTransparentClusterPosition(GetTransparentClusterSize()-1);
//		Int_t channel1 = getTransparentClusterPosition(GetTransparentClusterSize()-2);
//
//	}

	return clusterChannel.front();
}

/**
 * highest channel number of cluster...
 * goes to end of cluster and looks what is the channel no
 * @return highest channel number of cluster
 */
UInt_t TCluster::getHighestChannelNumber()
{
	if(clusterChannel.size()==0)
		return 0;
//	if(IsTransparentCluster())
//		return TMath::Max(getChannel(GetTransparentClusterSize()-2),getChannel(GetTransparentClusterSize()-1));
	return clusterChannel.back();
}


Float_t TCluster::getHighestSignal(bool cmnCorrected){
	Float_t highestSignal=0;
	UInt_t highestSignalClusterPos = -1;
	UInt_t highestSignalChannelPos;
	for(UInt_t clPos=0;clPos<checkClusterForSize();clPos++){
		Float_t signal = getSignal(clPos,cmnCorrected);
		if(signal>highestSignal){
			highestSignal=signal;
			highestSignalClusterPos=clPos;
		}
	}
	highestSignalChannelPos=getChannel(highestSignalClusterPos);
	if(maximumSignal!=highestSignal){
		cout<<"maximumSignal "<<maximumSignal<<" and highest signal "<<highestSignal<<" does not match... Something is wrong:";
		cout<<"highestSignalChannelPos: "<<highestSignalChannelPos<<"\thighesSignalClusterPos: "<<highestSignalClusterPos<<" "<<getHighestHitClusterPosition()<<endl;

	}

	return this->maximumSignal;
}



UInt_t TCluster::getClusterPosition(UInt_t channelNo){
	if(channelNo<this->getSmallestChannelNumber()&&channelNo>this->getHighestChannelNumber()){
		if(verbosity)cout<<"ChannelNo does not match..."<<endl;
		return 9999;
	}
	UInt_t clPos;
	for(clPos=0;clPos<checkClusterForSize()&&getChannel(clPos)!=channelNo;clPos++){}
	return clPos;
}

Float_t TCluster::getSignalOfChannel(UInt_t channel, bool cmnCorrected)
{
	if(channel<this->getSmallestChannelNumber()&&channel>this->getHighestSignalChannel()) return 0;
	UInt_t clPos = getClusterPosition(channel);
	if(clPos<getClusterSize()){//TODO
		Float_t signal = getSignal(clPos,cmnCorrected);
		//		if (signal>0)
		return signal;
	}
	//	cout << "TCluster::getSignalOfChannel: clPos>=getClusterSize()\treturn 0" << endl;
	return 0;
}


Float_t TCluster::getSignal(UInt_t clusterPos, bool cmnCorrected)
{
	if(clusterPos<checkClusterForSize()){
		Int_t adc = getAdcValue(clusterPos);
		Float_t pedMean = getPedestalMean(clusterPos,cmnCorrected);
		Float_t signal = (Float_t)adc - pedMean;
		if (cmnCorrected) signal -= getCMN();
		//		 if(!cmnCorrected && signal != getSignal(clusterPos))
		//			 cout<< "signal clalulated: "<<signal<<" other:"<<getSignal(clusterPos)<<endl;
		return signal;
	}
	else {
		if(verbosity)cout<<"clusterPos "<<clusterPos<<" bigger than clusterSize"<<checkClusterForSize()<<endl;
		return 0;
	}
}


/**
 */
Float_t TCluster::getSNR(UInt_t clusterPos, bool cmnCorrected)
{

	if(clusterPos<checkClusterForSize()){
		Float_t signal= getSignal(clusterPos,cmnCorrected);
		Float_t sigma = getPedestalSigma(clusterPos,cmnCorrected);
		return signal/sigma;
	}
	else return -1;

}
/**
 */
Float_t TCluster::getPedestalMean(UInt_t clusterPos, bool cmnCorrected)
{

	if(clusterPos<checkClusterForSize())
		if( cmnCorrected)
			return this->clusterPedMeanCMN.at(clusterPos);
		else
			return this->clusterPedMean.at(clusterPos);
	else return -1;

}


Float_t TCluster::getPedestalSigma(UInt_t clusterPos,bool cmnCorrected)
{
	if(clusterPos<checkClusterForSize())
		if( cmnCorrected)
			return this->clusterPedSigmaCMN.at(clusterPos);
		else
			return this->clusterPedSigma.at(clusterPos);
	else return -1;

}

Int_t TCluster::getAdcValue(UInt_t clusterPos)
{
	if(clusterPos<checkClusterForSize())
		return this->clusterADC.at(clusterPos);
	else return 0;
}

/**
 * function which is returning a variable which is similar to eta:
 *  S_L = signal which is in the channel  left from highest
 *  S_R = signal which is in the channel right from the highest
 *  S_H = signal from the highest signal
 *
 * @return S_L/(S_R+S_H)
 */
Float_t TCluster::getLeftEta(bool cmnCorrected){

	UInt_t clPosHighest = getHighestHitClusterPosition();
	UInt_t clPosLeft = clPosHighest-1;
	if(clPosLeft>=0){
		Float_t signalHighest = getSignal(clPosHighest,cmnCorrected);
		Float_t signalLeft = getSignal(clPosLeft,cmnCorrected);
		Float_t sumSignal = (signalLeft+signalHighest);
		if(sumSignal==0)
			return -1;
		return signalLeft/sumSignal;
	}
	return -1;
}

/**
 * function which is returning a variable which is similar to eta:
 *  S_L = signal which is in the channel  left from highest
 *  S_R = signal which is in the channel right from the highest
 *  S_H = signal from the highest signal
 *
 * @return S_R/(S_R+S_H)
 */
Float_t TCluster::getRightEta(bool cmnCorrected){

	UInt_t clPosHighest = getHighestHitClusterPosition();
	UInt_t clPosRight = clPosHighest+1;
	if(clPosRight<getClusterSize()){
		Float_t signalHighest = getSignal(clPosHighest,cmnCorrected);
		Float_t signalRight = getSignal(clPosRight,cmnCorrected);
		Float_t sumSignal = (signalRight+signalHighest);
		if(sumSignal==0)
			return -1;
		return signalRight/sumSignal;
	}
	return -1;
}

/**
 * @brief Calculation of reversed Eta of the cluster
 *	signalLeft/ (signalLeft+signalRight);
 * @return eta or -1 if not a valid Cluster
 */

Float_t TCluster::getReversedEta(Int_t &rightChannel,bool cmnCorrected){
	Int_t leftChannel;
	Float_t eta = getEta(leftChannel,cmnCorrected);
	if(eta>=0&&eta<=1){
		rightChannel = leftChannel+1;
		return 1-eta;
	}
	else{
		rightChannel = leftChannel;
	}
	return -1;
}

Float_t TCluster::getReversedEta(bool cmnCorrected){
	Int_t rightChannel=0;
	return getReversedEta(rightChannel,cmnCorrected);
}

Float_t TCluster::getEta(bool cmnCorrected){
	Int_t leftChannel=0;
	return getEta(leftChannel,cmnCorrected);
}


Float_t  TCluster::getEta(UInt_t clusPos,Int_t &leftChannel,bool cmnCorrected){

	UInt_t clPos2ndHighest = getHighestSignalNeighbourClusterPosition(clusPos,cmnCorrected);
	leftChannel = getChannel(clusPos);
	if(clPos2ndHighest>checkClusterForSize()){
		return -1;
	}
	UInt_t leftClPos = 0;
	UInt_t rightClPos = 0;
	if (clusPos < clPos2ndHighest) {
		leftClPos = clusPos;
		rightClPos = clPos2ndHighest;
	}
	else {
		leftClPos = clPos2ndHighest;
		rightClPos = clusPos;
	}
	leftChannel = this->getChannel(leftClPos);
	Float_t signalLeft  = getSignal(leftClPos,cmnCorrected);
	Float_t signalRight = getSignal(rightClPos,cmnCorrected);
	if(signalLeft<0||signalRight<0){
		if(signalRight>signalLeft)
			leftChannel = this->getChannel(rightClPos);
		return -1;
	}
	Float_t sumSignal = (signalLeft+signalRight);
	if(sumSignal==0){
		return -1;
	}

	Float_t eta = signalRight/sumSignal;
	if(verbosity&&(eta<0||eta>1))
		this->Print(1);
	return eta;
}
/**
 * @brief Calculation of Eta of the cluster
 * 			signalRight / (signalLeft+signalRight);
 * @return eta, -2 if clusterSize<1, -1 if invalid cluster
 */
Float_t TCluster::getEta(Int_t &leftChannel,bool cmnCorrected)
{
	leftChannel = -1;
	if (checkClusterForSize() < 1) return -2;
	if(checkClusterForSize()==1){
		leftChannel= getChannel(0);
		return -1;
	}
	if (IsTransparentCluster()){
		if( GetTransparentClusterSize() < 2){
			leftChannel = getChannel(getHighestHitClusterPosition());
			return -1;
		}
	}
	return getEta(getHighestHitClusterPosition(),leftChannel,cmnCorrected);
//	UInt_t clPosHighest = getHighestHitClusterPosition();
//	UInt_t clPos2ndHighest = getHighestSignalNeighbourClusterPosition(clPosHighest);
//	UInt_t leftClPos = 0;
//	UInt_t rightClPos = 0;
//	if (clPosHighest < clPos2ndHighest) {
//		leftClPos = clPosHighest;
//		rightClPos = clPos2ndHighest;
//	}
//	else {
//		leftClPos = clPos2ndHighest;
//		rightClPos = clPosHighest;
//	}
//	leftChannel = this->getChannel(leftClPos);
//	Float_t signalLeft  = getSignal(leftClPos,cmnCorrected);
//	Float_t signalRight = getSignal(rightClPos,cmnCorrected);
//	Float_t sumSignal = (signalLeft+signalRight);
//	if(sumSignal==0)
//		return -1;
//	Float_t eta = signalRight/sumSignal;
//	if(verbosity&&(eta<0||eta>1))
//		this->Print(1);
//	return eta;
}

Float_t TCluster::getEtaPostion(bool cmnCorrected){
	Int_t leftChannel;
	Float_t eta = getEta(leftChannel,cmnCorrected);
	if(leftChannel==-1)return INVALID_POSITION;
	if(eta<0)
		return leftChannel;
	return eta+leftChannel;
}

Float_t TCluster::getPositionCorEta(bool cmnCorrected,TH1F* histo){
	Int_t leftChannel;
	Float_t eta = getEta(leftChannel,cmnCorrected);
	if(eta<=0){
		return getEtaPostion(cmnCorrected);
	}
	if(eta>1){
		return getEtaPostion(cmnCorrected);
	}
	Float_t corEta= getValueOfHisto(eta,histo);
	if(verbosity)	cout<<leftChannel<<" + "<<corEta<<" = "<<leftChannel+corEta;
	if(leftChannel==-1)
		return INVALID_POSITION;
	return leftChannel+corEta;
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

UInt_t TCluster::getFirstHitChannel(){
	UInt_t cl = 0;
	UInt_t size =  this->getClusterSize();
	if (size==0) return -1;
	while (!isHit(cl) && cl < size)
		cl++;
	return getChannel(cl);
}

UInt_t TCluster::getLastHitChannel(){
    if (this->getClusterSize()) return -1;
	UInt_t cl = this->getClusterSize()-1;

	while ( !isHit(cl) && cl >= 0)
		cl--;
	return getChannel(cl);
}

UInt_t TCluster::getHighestSignalNeighbourChannel(UInt_t channelNo,bool cmnCorrected)
{
	return getChannel(getHighestSignalNeighbourClusterPosition(getClusterPosition(channelNo),cmnCorrected));
}

/**
 * gets the clusterposition of the highest signal neigboured to a given clusterposition
 * @param clPos position where should be looked next to
 * @return clusterPosition in cluster
 */
UInt_t TCluster::getHighestSignalNeighbourClusterPosition(UInt_t clPos,bool cmnCorrected,bool bNegativeSignals)
{
    if (false&&bNegativeSignals)
        cout <<"TCluster::getHighestSignalNeighbourClusterPosition"<<clPos<<"/"<<cmnCorrected<<
        "/clSize: "<<checkClusterForSize()<<endl;
	if (clPos>=checkClusterForSize() || clPos<0 || checkClusterForSize()<2) return 9999;
	if(checkClusterForSize()==2){
		if(clPos==1) return clPos-1;
		else if(clPos==0) return clPos+1;
		else return 9999;
	}
	if(IsTransparentCluster()){

		bool valid1 = IsValidTransparentClusterPosition(clPos-1);
		bool valid2 = IsValidTransparentClusterPosition(clPos+1);
		//cout<<"Valid Hits: "<<clPos-1<<": "<<valid1<<"/ "<<clPos+1<<": "<<valid2<<endl;
		if(!valid1&&!valid2)
			return 9999;
		else if (!valid1)
			return clPos+1;
		else if (!valid2)
			return clPos-1;
	}
	Float_t signalLeft = getSignal(clPos-1,cmnCorrected);
	Float_t signalRight=getSignal(clPos+1,cmnCorrected);
	if (false&&bNegativeSignals)
	    cout<<"clPos-1: "<<signalLeft<<"\tclPos+1:"<<signalRight<<endl;
	if (signalLeft < signalRight){
		if (signalRight>0|| bNegativeSignals)
			return clPos+1;
	}
	else{
		if (signalLeft > 0||bNegativeSignals) return clPos-1;
	}
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
	if(histo==0){
		cerr<<"Histo pointer = 0!"<<endl;
		return -999;
	}

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


bool TCluster::hasInvalidReadout(){
//	cout<<"\t\tTCluster::hasInvalidReadout\t "<<flush;
	bool invalidReadout = false;
	Float_t minInvalidSignal = TPlaneProperties::GetMinInvalidSignal(det);
	UInt_t cl;
	for(cl=0;cl<checkClusterForSize()&&!invalidReadout;cl++){
//		cout<<cl<<" / "<<checkClusterForSize()<<flush;
		invalidReadout = getSignal(cl)<minInvalidSignal||invalidReadout;
	}
//	cout<<" last "<<cl<<"__\t__"<<flush;
//	if(invalidReadout)
//		cout<<"InvalidReadout "<<cl<<flush;
//	cout<<endl;
	return invalidReadout;
}

void TCluster::Print(UInt_t level){
	if(level>10)cout<<"\n";
	cout<<Intent(level%10)<<"Cluster of Event "<<flush;
	cout<<eventNumber<<", det "<<(int)det<<" with "<<size()<<"/"<<checkClusterForSize()<<" Ch, CMN "<<cmNoise;
	if(IsTransparentCluster())
		cout<<" transCluster @ "<<isTransparentCluster<<" size: "<<GetTransparentClusterSize();
//	cout<<": "<<getCharge(false,true)<<"/"<<getCharge(true,true)<<"\n\t"<<flush;
	for(UInt_t cl=0;cl<checkClusterForSize();cl++){
		if(this->isSeed(cl))
			cout<<"\t{";
		else if(this->isHit(cl))
			cout<<"\t(";
		else
			cout<<"\t[";
		cout<<this->getChannel(cl)<<"|"<<this->getAdcValue(cl)<<"|"<<this->getSignal(cl,false)<<","<<this->getSignal(cl,true)<<"|"<<this->getSNR(cl);
		if(this->isSeed(cl))
			cout<<"}"<<flush;
		else if(this->isHit(cl))
			cout<<")"<<flush;
		else
			cout<<"]"<<flush;
	}
	cout<<"\t||"<<flush;
	cout<<this->getHighestSignalChannel()<<" "<<flush<<this->getHighest2Centroid()<<" "<<this->getChargeWeightedMean(true)<<" "<<this->getEtaPostion(false)<<" "<<this->getPositionCorEta(false);
	cout<<" "<<this->getEta()<<" "<<this->getEta(true);
	cout<<endl;
}

Float_t TCluster::getPositiveCharge(UInt_t nClusterEntries, bool cmnCorrected,
        bool useSmallSignals) {
    if(nClusterEntries==0)return 0;
    Int_t clusPosStart = getClusterPosition(getHighestSignalChannel());
    direction_t dir = getSignal(clusPosStart-1,cmnCorrected)>getSignal(clusPosStart+1,cmnCorrected)?left:right;
    Float_t clusterCharge = getChargeStartingAt(nClusterEntries,clusPosStart,dir,true,cmnCorrected,useSmallSignals);
    return clusterCharge;
}
/**
 *
 * @param charge
 * @param pos
 * @param cmnCorrected
 * @return
 */
bool TCluster::hasNegativeCharge(Float_t& charge, Int_t& pos, bool cmnCorrected, bool verb) {
    Float_t oldCharge = 0;
    Float_t currentCharge;
    if (verb) cout<<"[TCluster::hasNegativeCharge]: "<<endl;
    bool hasNegCharge = false;
    if (verb) this->Print();
    pos = 0;
    Float_t negCharge = 0;
    Float_t smallCharge = 1e9;

    for(UInt_t clusterSize = 1; clusterSize<= this->GetTransparentClusterSize();clusterSize++){
        currentCharge = this->getCharge(clusterSize,cmnCorrected);
        charge = currentCharge - oldCharge;
        Int_t clusterNo = clusterSize -1;
        Int_t clPos = this->getTransparentClusterPosition(clusterNo);
        Int_t startChannel = Int_t(isTransparentCluster+.5);
        Int_t clStart = this->getClusterPosition(startChannel);
        Int_t orig_dir;
        Float_t delta = (Float_t)startChannel-isTransparentCluster;
        if (delta<0)
            orig_dir = -1;
        else
            orig_dir = 1;
        Int_t dir = orig_dir;
        if (clusterNo%2==0) dir *= -1;
        Int_t dif = (Int_t(clusterNo)+1)/2;
        Int_t clusPos = clStart + dir *dif;
        if (clusPos != clPos){
            cout<<"\n\nclStart: "<<clStart<<" origDir: "<<orig_dir<<" dir: "<<dir<<" dif: "<<dif<<" clusPos:"<<clusPos<<"/"<<clPos<<endl;
            this->Print();
            char t;
            cin>>t;
        }

        Float_t signal = this->getSignal(clPos,cmnCorrected);
        if (charge != signal && !verb){
            cout<<"\n\n"<<clusterSize<<"/"<<clPos<<"/"<<clusterNo<<": Something is wrong: "<<charge<<"/"<<signal<<endl;
            cout<<"clusNo: "<<clusterNo<<"startCHannel: "<<startChannel<<" clStart: "<<clStart
                    <<" Delta: "<<(startChannel-isTransparentCluster)
                <<" origDir: "<<orig_dir<<" dir: "<<dir<<" dif: "<<dif<<" clusPos:"<<clusPos<<"/"<<clPos<<endl;
            return hasNegativeCharge(charge,pos,cmnCorrected,true);
            verb = true;
        }
        if (charge < smallCharge){
            smallCharge = charge;
        }
        if (verb) cout<<"\t"<<clusterSize<<"-"<<clPos<<TString::Format("%+7.1f - %+7.1f - %+7.1f | %7.1f | %7.1f",
                currentCharge,oldCharge,charge,smallCharge, signal)<<"\n";
        if (charge < 0 and !hasNegCharge){
            if (clusterSize==1 && clusterSize+2<=this->size()){
                pos = 0;
                hasNegCharge = true;
                negCharge = charge;
                Int_t pos1 = this->getTransparentClusterPosition(clusterSize);
                Int_t pos2 = this->getTransparentClusterPosition(clusterSize+1);
                Float_t sig1 = this->getSignal(pos1,cmnCorrected);
                Float_t sig2 = this->getSignal(pos2,cmnCorrected);
                bool change = false;
                if (sig1<sig2){
                    if (sig1<charge){
                        negCharge = sig1;
                        pos = 1;
                        change = true;
                    }
                }
                else{
                    if (sig2<charge){
                        negCharge = sig2;
                        pos = 2;
                        change=true;
                    }
                }
                if (verb||change) 
                    cout <<clusterSize<< " found negative charge at "<< pos<<": "<<signal<<"/"<<charge<<"/"
                        <<"/"<<sig1<<"/"<<sig2<<endl;
            }
            else if( clusterSize+1<this->size() &&
                this->getCharge(clusterSize+1,cmnCorrected)-currentCharge > charge){
                pos = clusterSize;
                if (verb) cout <<clusterSize<< " found negative charge at "<< pos<<": "<<charge<<endl;
                hasNegCharge = true;
                negCharge = charge;
            }else{
                pos = clusterSize+1;
                charge =  this->getCharge(clusterSize+1,cmnCorrected)-currentCharge;
                if (verb) cout << clusterSize<<" found negative charge at "<< pos<<": "<<charge<<endl;
                hasNegCharge =  true;
                negCharge = charge;
            }
        }
        oldCharge = currentCharge;
    }
    if (false && pos !=3 && hasNegCharge){
        cout<<"press a key"<<endl;
        char t;
        cin>>t;
    }
    if (hasNegCharge)
        charge = negCharge;
    else
        charge = smallCharge;
    if(hasNegCharge && verb)
        cout<<"found negative charge at "<< pos<<": "<<negCharge<<endl;
    else if (verb)
        cout<<"no negative charge at "<< pos<<": "<<negCharge<<endl;

    return hasNegCharge;
}


