/*
 * TTracking.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: bachmair
 */

#include "../include/TTracking.hh"

TTracking::TTracking(std::string pathName, std::string alignmentName,std::string etaDistributionPath, UInt_t runNumber):TADCEventReader(pathName,runNumber){
	alignmentFile=NULL;
	cout<<"new TTracking: \n\tpathName:"<<pathName<<"\n\ta;ignmentName: "<<alignmentName<<"\n\tetaDistPath: "<<etaDistributionPath<<"\n\tRunNumber: "<<runNumber<<endl;
	setAlignment(alignmentName);
	this->setEtaDistributionPath(etaDistributionPath);
	if(myAlignment==0){
		cerr<<"Not able to Read Alignment. EXIT!"<<endl;
		exit(-1);
	}

	if(myAlignment!=NULL)
		myTrack=new TTrack(myAlignment);
	else
		myTrack=NULL;
	if(myTrack!=NULL)
		for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
		{
      cout<<"Set EtaIntegral of detector "<<det<<flush;
      TH1F* etaIntegral=(TH1F*)this->getEtaIntegral(det);
			myTrack->setEtaIntegral(det,etaIntegral);
			cout<<" successful"<<endl;
		}
	cout<<"DONE WITH LOADING TTRACKING"<<endl;
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
	return myTrack->predictPosition(subjectPlane,vecRefPlanes,TCluster::corEta,bPrint);
}

Float_t TTracking::getXPosition(UInt_t plane)
{
	if(myTrack==0)
			return 0;
		return myTrack->getXPosition(plane);
}

Float_t TTracking::getYPosition(UInt_t plane)
{
	if(myTrack==0)
			return 0;
		return myTrack->getYPosition(plane);
}

Float_t TTracking::getZPosition(UInt_t plane)
{
	if(myTrack==0)
			return 0;
		return myTrack->getZPosition(plane);
}

Float_t TTracking::getMeasured(TPlaneProperties::enumCoordinate cor, UInt_t plane,TCluster::calculationMode_t mode)
{
	if (myTrack == 0)
		return 0;
	return myTrack->getMeasured(cor,plane,mode);
}

bool TTracking::setAlignment(std::string alignmentName){
	if(this->alignmentFile!=NULL) alignmentFile->Delete();
	alignmentFile=NULL;
	cout<<"load AlignmentFile: \""<<alignmentName<<"\""<<endl;
	alignmentFile = new TFile(alignmentName.c_str());
	alignmentFile->GetObject("alignment",myAlignment);
	if(myAlignment==NULL){
		cerr<<"COULD NOT READ THE ALIGNMENT FILE..."<<endl;
		return false;
	}
	else{
		cout<<"Read the Alignment of AlignmentFile...\n"<<endl;
		myAlignment->PrintResults();
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
    cout<<"eventNumber > entries: "<<eventNumber<<" "<<GetEntries()<<endl;
    eventNumber=GetEntries()-1;
    if(eventNumber<0)eventNumber=0;
    cout<<"new eventNumber ==> "<<eventNumber;
  }
  if(myTrack!=NULL){
    if(verbosity>3){cout<<"myTrack!=0"<<endl;};
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
	return myTrack->getPositionInLabFrameStripDetector(plane,xCluster,yPred,mode,histo);
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
Float_t  TTracking::getPosition(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionInLabFrame(cor,plane,mode);
}
Float_t TTracking::getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred){
	if(myTrack==0)
		return 0;
	return myTrack->getPositionInDetSystem(det,xPred,yPred);
}
