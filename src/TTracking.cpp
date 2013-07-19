/*
 * TTracking.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: bachmair
 */

#include "../include/TTracking.hh"

TTracking::TTracking(std::string pathName, std::string alignmentName,std::string etaDistributionPath, TSettings*settings):TADCEventReader(pathName,settings){
	this->settings= settings;
	alignmentFile=NULL;
	verbosity = settings->getVerbosity();
	if(verbosity) cout<<"new TTracking: \n\tpathName:"<<pathName<<"\n\ta;ignmentName: "<<alignmentName<<"\n\tetaDistPath: "<<etaDistributionPath<<"\n\tRunNumber: "<<settings->getRunNumber()<<endl;
	setAlignment(alignmentName);
	this->setEtaDistributionPath(etaDistributionPath);
	if(myAlignment==0){
		cerr<<"Not able to Read Alignment. EXIT!"<<endl;
		exit(-1);
	}

	if(myAlignment!=NULL)
		myTrack=new TTrack(myAlignment,this->settings);
	else
		myTrack=NULL;
	if(myTrack!=NULL)
		for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		{
			if(verbosity) cout<<"Set EtaIntegral of detector "<<det<<flush;
			TH1F* etaIntegral=(TH1F*)this->getEtaIntegral(det);
			myTrack->setEtaIntegral(det,etaIntegral);
			if(verbosity) cout<<" successful"<<endl;
		}
	if(verbosity) cout<<"DONE WITH LOADING TTRACKING"<<endl;
}

TTracking::~TTracking() {
	delete myTrack;
	delete myAlignment;
	delete alignmentFile;
}

TPositionPrediction *TTracking::predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes, bool bPrint)
{
	if(myTrack==0)
		return 0;
	if (this->getEvent() ==0)
		cerr<<"pEvent pointer is null...."<<endl;
	myTrack->setEvent(this->getEvent());
	return myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,bPrint);
}

Float_t TTracking::getXPosition(UInt_t plane)
{
	if(myTrack==0)
		return 0;
	return myTrack->getXPositionMetric(plane);
}

Float_t TTracking::getYPosition(UInt_t plane)
{
	if(myTrack==0)
		return 0;
	return myTrack->getYPositionMetric(plane);
}

Float_t TTracking::getZPosition(UInt_t plane)
{
	if(myTrack==0)
		return 0;
	return myTrack->getZPosition(plane);
}

Float_t TTracking::getMeasuredPositionMetricSpace(TPlaneProperties::enumCoordinate cor, UInt_t plane,TCluster::calculationMode_t mode)
{
	if (myTrack == 0)
		return 0;
	return myTrack->getMeasuredClusterPositionMetricSpace(cor,plane,mode);
}

Float_t TTracking::getYPositionInDetSystem(UInt_t det, Float_t xPredicted,
		Float_t yPredicted) {
//	UInt_t subjectPlane = TPlaneProperties::getDiamondPlane();
//	UInt_t subjectDetector = TPlaneProperties::getDetDiamond();
//
//	predictedPosition = eventReader->predictPosition(subjectPlane,vecSilPlanes);
	UInt_t subjectPlane = TPlaneProperties::getPlaneNumber(det);
	Double_t xOffset = myAlignment->GetXOffset(subjectPlane);
	Double_t PhiOffset = myAlignment->GetPhiXOffset(subjectPlane);

	Float_t yDet = yPredicted*TMath::Cos(PhiOffset)+(xPredicted-xOffset)*TMath::Sin(PhiOffset);
	yDet -= settings->get3DYOffset();
	return yDet;
}

bool TTracking::setAlignment(std::string alignmentName){
	if(this->alignmentFile!=NULL) alignmentFile->Delete();
	alignmentFile=NULL;
	if(verbosity) cout <<"load AlignmentFile: \""<<alignmentName<<"\""<<endl;
	alignmentFile = new TFile(alignmentName.c_str());
	alignmentFile->GetObject("alignment",myAlignment);
	if(myAlignment==NULL){
		cerr<<"COULD NOT READ THE ALIGNMENT FILE..."<<endl;
		return false;
	}
	else{
		if(verbosity) cout<<"Read the Alignment of AlignmentFile...\n"<<endl;
		if(verbosity) myAlignment->PrintResults();
		return true;
	}
	return true;
}

bool TTracking::LoadEvent(UInt_t eventNumber){
	if(verbosity>3){
		//    cout<<"Load Event: "<<eventNumber<<flush;
		//    cout<<" "<<GetEntries()<<endl;
	}

	if(eventNumber>GetEntries()){
		if(verbosity) cout<<"eventNumber > entries: "<<eventNumber<<" "<<GetEntries()<<endl;
		eventNumber=GetEntries()-1;
		if(eventNumber<0)eventNumber=0;
		if(verbosity) cout<<"new eventNumber ==> "<<eventNumber;
	}
	if(myTrack!=NULL){
		if(verbosity>6){cout<<"myTrack!=0"<<endl;};
		bool retVal=TADCEventReader::LoadEvent(eventNumber);
		if(retVal)
			myTrack->setEvent(this->getEvent());
		return retVal;
	}
	else{
		cout<<"MYTRACK==0!!!!!! WHATS WRONG???????"<<endl;
		char t;cin>>t;
		exit(-1);
	}

	return false;
}

Float_t  TTracking::getStripXPositionOfCluster(UInt_t plane,TCluster xCluster, Float_t yPred,TCluster::calculationMode_t mode,TH1F* histo){
	if(myTrack==0)
		return 0;
	if (histo==0)
		histo=getEtaIntegral(plane*2);
	return myTrack->getXPositionInLabFrameStripDetector(plane,xCluster,yPred,mode,histo);
}
Float_t  TTracking::getStripXPosition(UInt_t plane,Float_t yPred,TCluster::calculationMode_t mode){
	if(myTrack==0)
		return 0;
	TH1F* histo=getEtaIntegral(plane*2);
	return myTrack->getStripXPosition(plane,yPred,mode,histo);

}
Float_t  TTracking::getPositionOfCluster(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster xCluster,TCluster yCluster, TCluster::calculationMode_t mode, TH1F* histo){
	if(myTrack==0)
		return 0;
	return myTrack->getPostionInLabFrame(cor,plane,xCluster,yCluster,mode,histo);
}


Float_t TTracking::getPositionOfCluster(UInt_t det, TCluster cluster, Float_t predictedPerpPosition, TCluster::calculationMode_t mode, TH1F* histo){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionOfCluster(det,cluster,predictedPerpPosition,mode,histo);
}

/** from the current activated event positionInLabFrame
 *
 * @param cor
 * @param plane
 * @param mode
 * @return
 */
Float_t  TTracking::getPosition(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionInLabFrame(cor,plane,mode);
}

/**  from pediction calculated hit Position of a detector in detector system, metric space
 *
 * @param det
 * @param xPred
 * @param yPred
 * @return
 */
Float_t TTracking::getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionInDetSystem(det,xPred,yPred);
}
