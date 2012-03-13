//
//  TTrack.cpp
//  Diamond Analysis
//
//  Created by Lukas Bäni on 06.12.11.
//  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
//

#include "../include/TTrack.hh"

/**
 * @brief cornsturcter of TTrack class
 * @param alignment pointer on current alignment
 */
TTrack::TTrack(TDetectorAlignment *alignment) {
	// TODO: check if event has the one & only one cluster flag
	verbosity=0;
	this->alignment = alignment;
	event=NULL;

	linFitX = new TLinearFitter(1,"pol1","D");
	linFitY= new TLinearFitter(1,"pol1","D");
	linFitX->StoreData(true);
	linFitY->StoreData(true);

	linFitX->SetFormula("pol1");
	linFitY->SetFormula("pol1");
}

TTrack::~TTrack() {
	delete linFitX;
	delete linFitY;
}


/**
 * @brief gives  number of cluster in detector
 * @param det detector number
 * @return number of cluster in detector det
 */
UInt_t TTrack::getNClusters(int det) {
	if(event==NULL)return N_INVALID;
	if (det%2 == 0) {
		int plane = det / 2;
		return this->event->getPlane(plane).getNXClusters();
	}
	else {
		int plane = (det-1) / 2;
		return this->event->getPlane(plane).getNYClusters();
	}
}


/**
 * calculate Position in absoltue space of two cluster Hits in X/Y with calculation mode mode
 *
 * @param cor coordinate in which the calculation is made
 * @param plane which plane no is used
 * @param xCluster
 * @param yCluster
 * @param mode
 * @return value of calculated position
 */
Float_t TTrack::getPositionOfCluster(TPlane::enumCoordinate cor,UInt_t plane,TCluster xCluster,TCluster yCluster, TCluster::calculationMode_t mode,TH1F* histo){
	if(xCluster.size()<=0||yCluster.size()<=0)
		return N_INVALID;
	Float_t xOffset = this->getXOffset(plane);
	Float_t yOffset = this->getYOffset(plane);
	Float_t phiXOffset = this->getPhiXOffset(plane);
	Float_t phiYOffset = this->getPhiYOffset(plane);
	Float_t xMeasured = xCluster.getPosition(mode,getEtaIntergral(plane*2));
	Float_t yMeasured = yCluster.getPosition(mode,getEtaIntergral(plane*2+1));
	// apply offsets
	Float_t xPosition = (xMeasured) * TMath::Cos(phiXOffset) + (yMeasured) * TMath::Sin(phiYOffset);
	Float_t yPosition = (xMeasured) * TMath::Sin(-phiXOffset) + (yMeasured) * TMath::Cos(phiYOffset);
	xPosition += xOffset;
	yPosition += yOffset;
	switch(cor){
		case TPlane::X_COR: return xPosition;break;
		case TPlane::Y_COR: return yPosition;break;
		case TPlane::Z_COR: return getZPosition(plane);break;
		default: return N_INVALID;
	}
}

/**
 * calculate Position in lab system of a cluster
 *
 * @param det detector in which the cluster is
 * @param cluster
 * @param predictedPerpPosition predicted position of the track in the perpendicular det of the same plane (in lab system)
 * @param mode
 * @return value of calculated position
 */
Float_t TTrack::getPositionOfCluster(UInt_t det, TCluster cluster, Float_t predictedPerpPosition, TCluster::calculationMode_t mode,TH1F* histo) {
	if (cluster.size()<=0) return N_INVALID;
	UInt_t plane = det/2;
	if(histo==0)
		histo = (TH1F*)getEtaIntergral(plane*2+1);
	Float_t measuredPos = cluster.getPosition(mode,histo);
	Float_t xOffset = this->getXOffset(plane);
	Float_t yOffset = this->getYOffset(plane);
	Float_t phiXOffset = this->getPhiXOffset(plane);
	Float_t phiYOffset = this->getPhiYOffset(plane);
	Float_t xPosition = (measuredPos) / TMath::Cos(phiXOffset) + (predictedPerpPosition) * TMath::Tan(phiXOffset);
	Float_t yPosition = (measuredPos) / TMath::Cos(phiYOffset) - (predictedPerpPosition) * TMath::Tan(phiYOffset);
	xPosition += xOffset;
	yPosition += yOffset;
	if (det%2 == 0) {
		return xPosition;
	}
	else {
		return yPosition;
	}
}

/**
 * Calculation of Hitposition of event using the measured Offsets
 * in this transformation it transform the measured hitPosition in the
 * plane space into the space of the first plane
 * @param cor enum of Coordinate from which you want to have the result
 * @param plane number of plane for which you are calculating the hitposition
 * @return calculated hit position
 */
Float_t TTrack::getStripXPosition(UInt_t plane,Float_t yPred,TCluster::calculationMode_t mode,TH1F* histo){
	if(event==NULL)return N_INVALID;
	if(event->getNXClusters(plane)!=1)
		return N_INVALID;
	// get offsets
	TCluster xCluster = event->getCluster(plane,TPlane::X_COR,0);
	if(xCluster.size()<=0) {
		cerr<<"This cluster is too small!!!!"<<endl;
		xCluster.Print();
	}
	if(histo==NULL)
		histo = getEtaIntergral(plane*2);
	return getStripXPositionOfCluster(plane,xCluster,yPred,mode,histo);
}

/**
 * todo anpassen!!!
 * Calculation of Hitposition of event using the measured Offsets
 * in this transformation it transform the measured hitPosition in the
 * plane space into the space of the first plane
 * @param cor enum of Coordinate from which you want to have the result
 * @param plane number of plane for which you are calculating the hitposition
 * @return calculated hit position
 */
Float_t TTrack::getStripXPositionOfCluster(UInt_t plane,TCluster xCluster, Float_t yPred,TCluster::calculationMode_t mode,TH1F* histo){
	if(xCluster.size()<=0)
		return N_INVALID;
	// get offsets
	Float_t xOffset = this->getXOffset(plane);
	Float_t phiXOffset = this->getPhiXOffset(plane);
	if(histo==NULL)
		histo = getEtaIntergral(plane*2);

	Float_t xMeasured = xCluster.getPosition(mode,histo);//-xOffset;

	// apply offsets
	Float_t xPosition = (xMeasured) / TMath::Cos(phiXOffset) + (yPred) * TMath::Tan(phiXOffset);
	xPosition += xOffset;
	if(verbosity)cout<<"Xpos of"<<plane<<": with "<<phiXOffset<<" and  "<<xOffset<<" measured: "<<xMeasured<<", "<<xPosition<<"/"<<yPred<<endl;

	return xPosition;
}
/**
 * gives the xPosition of "event" in the telescope space
 * this function uses getPosition(X_COR,plane)
 * @param plane planeNumber to use correct offsets
 * @return calculated xPosition
 */
Float_t TTrack::getXPosition(UInt_t plane,TCluster::calculationMode_t mode,TH1F* histo) {
	if(histo==0)
		histo = getEtaIntergral(plane*2);
	return getPosition(TPlane::X_COR,plane,mode,histo);
}

/**
 * gives the yPosition of "event" in the telescope space
 * this function uses getPosition(Y_COR,plane)
 * @param plane planeNumber to use correct offsets
 * @return calculated yPosition
 */
Float_t TTrack::getYPosition(UInt_t plane,TCluster::calculationMode_t mode,TH1F* histo) {
	if(histo==0)
			histo = getEtaIntergral(plane*2);
	return getPosition(TPlane::Y_COR,plane,mode,histo);
}

Float_t TTrack::getZPosition(UInt_t plane,TCluster::calculationMode_t mode){
	return alignment->GetZOffset(plane);
}



TPositionPrediction* TTrack::predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,TCluster::calculationMode_t mode,bool bPrint)
{
	linFitX->ClearPoints();
	linFitY->ClearPoints();
	if(event==NULL){
		cerr<<"TTrack:predictPosition no ReferencePlanes are defined...event =NULL"<<endl;
		TPositionPrediction* prediction=0;
		return prediction;
	}
	if(vecRefPlanes.size()==0){
		cerr<<"TTrack:predictPosition no ReferencePlanes are defined...vecRefSize=0"<<endl;
		TPositionPrediction *prediction=0;
		return prediction;
	}
	if(vecRefPlanes.size()==1){
		if(verbosity>3)	cout<<"TTrack::predictPosition with 1 refPlane"<<endl;
		TPositionPrediction *prediction=new TPositionPrediction(getXPosition(vecRefPlanes.at(0),mode), 0.,0.,getYPosition(vecRefPlanes.at(0),mode),0.,0.);
		return prediction;
	}
	vector<Double_t> zPosVec;//todo add xsigma ysigma
	if(bPrint)cout<<"Prediction of Track in Plane "<<subjectPlane<<" with "<<vecRefPlanes.size()<<" Planes:"<<endl;
	for(UInt_t pl=0;pl<vecRefPlanes.size();pl++){
		UInt_t plane=vecRefPlanes.at(pl);
		zPosVec.clear();
		zPosVec.push_back(alignment->GetZOffset(plane));
		linFitX->AddPoint(&zPosVec.at(0),(Double_t)getXPosition(plane,mode),this->alignment->getXResolution(plane));//todo anpassen des SIGMA
		linFitY->AddPoint(&zPosVec.at(0),(Double_t)getYPosition(plane,mode),this->alignment->getYResolution(plane));//todo anpassen des sigma 0.001
		if(verbosity>3||bPrint)	cout<<"\tAdd in Plane "<<plane<<"  "<<getXPosition(plane,mode)<<"+/-"<<alignment->getXResolution(plane)<<"/"<<getYPosition(plane,mode)<<"+/-"<<alignment->getYResolution(plane)<<"/"<<getZPosition(plane)<<endl;
	}
	linFitX->Eval();
	linFitY->Eval();
	linFitX->Chisquare();
	linFitY->Chisquare();
	Float_t zPos = alignment->GetZOffset(subjectPlane);
	Float_t zSigma = 0; //todo
	Float_t mx = linFitX->GetParameter(1);
	Float_t sigma_mx = linFitX->GetParError(1);
	Float_t bx = linFitX->GetParameter(0);
	Float_t sigma_bx = linFitX->GetParError(0);
	Float_t my = linFitY->GetParameter(1);
	Float_t sigma_my = linFitY->GetParError(1);
	Float_t by = linFitY->GetParameter(0);
	Float_t sigma_by = linFitY->GetParError(0);
	Float_t xChi2 = linFitX->GetChisquare()/linFitX->GetNumberFreeParameters();
	Float_t yChi2 = linFitY->GetChisquare()/linFitY->GetNumberFreeParameters();
	Float_t xPos = mx*zPos+bx;
	Float_t yPos = my*zPos+by;
	Float_t xSigma = (zPos*sigma_mx)*(zPos*sigma_mx)+(mx*zSigma)*(mx*zSigma)+(sigma_bx*sigma_bx);
	xSigma = TMath::Sqrt(xSigma);
	Float_t ySigma = (zPos*sigma_my)*(zPos*sigma_my)+(my*zSigma)*(my*zSigma)+(sigma_by*sigma_by);
	ySigma = TMath::Sqrt(ySigma);
	TPositionPrediction* prediction=new TPositionPrediction(xPos,xSigma,xChi2,yPos,ySigma,yChi2);
	if(verbosity>3||bPrint)	cout<<"\n  Predition of Plane "<<subjectPlane<<" with "<<vecRefPlanes.size()<<"Planes: ZPosition: "<<zPos<<endl;
	if(verbosity>3||bPrint)	cout<<"\tX: "<<xPos<<" +/- "<<xSigma<<"   with a Chi^2 of "<<xChi2<<"  "<<linFitX->GetNpoints()<<endl;
	if(verbosity>3||bPrint)	cout<<"\tY: "<<yPos<<" +/- "<<ySigma<<"   with a Chi^2 of "<<yChi2<<"  "<<linFitY->GetNpoints()<<"\n"<<endl;
	return prediction;
}

vector<Float_t> TTrack::getSiXPositions() {
	vector<Float_t> xPositions;
	if(event==NULL)return xPositions;
	for (int plane = 0; plane < 4; plane++) {
		xPositions.push_back(this->getXPosition(plane));
	}
	return xPositions;
}

vector<Float_t> TTrack::getSiYPositions() {

	vector<Float_t> yPositions;
	if(event==NULL)return yPositions;
	for (int plane; plane < 4; plane++) {
		yPositions.push_back(this->getYPosition(plane));
	}
	return yPositions;
}

void TTrack::setDetectorAlignment(TDetectorAlignment *alignment)
{
	if(alignment!=this->alignment)
		delete this->alignment;
	this->alignment=alignment;
}

Float_t TTrack::getXMeasured(UInt_t plane,TCluster::calculationMode_t mode)
{
	return getMeasured(TPlane::X_COR,plane,mode);
}

Float_t TTrack::getYMeasured(UInt_t plane,TCluster::calculationMode_t mode)
{
	return getMeasured(TPlane::Y_COR,plane,mode);
}
/**
 * Calculation of Hitposition of event using the measured Offsets
 * in this transformation it transform the measured hitPosition in the
 * plane space into the space of the first plane
 * @param cor enum of Coordinate from which you want to have the result
 * @param plane number of plane for which you are calculating the hitposition
 * @return calculated hit position
 */
Float_t TTrack::getPosition(TPlane::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode,TH1F* histo){
	if(event==NULL)return N_INVALID;
	if(event->getNXClusters(plane)!=1||event->getNYClusters(plane)!=1)
		return N_INVALID;
	// get offsets
	if(histo==0)
			histo = getEtaIntergral(plane*2);
	TCluster xCluster,yCluster;
	xCluster=event->getPlane(plane).getXCluster(0);
	yCluster=event->getPlane(plane).getYCluster(0);
	return getPositionOfCluster(cor,plane,xCluster,yCluster,mode,histo);
}


Float_t TTrack::getMeasured(TPlane::enumCoordinate cor, UInt_t plane,TCluster::calculationMode_t mode)
{
	if(event==NULL)return N_INVALID;
	if(cor==TPlane::XY_COR&&(event->getNXClusters(plane)!=1||event->getNYClusters(plane)!=1))
	return N_INVALID;
	if(cor==TPlane::X_COR&&event->getNXClusters(plane)==1)
		return event->getPlane(plane).getXPosition(0,mode);
	if(cor==TPlane::Y_COR&&event->getNYClusters(plane)==1)
		return event->getPlane(plane).getYPosition(0,mode);
	return N_INVALID;
//// get offsets
//	switch(cor){
//	case TPlane::X_COR:break;
//	case TPlane::Y_COR:break;
//	default: return N_INVALID;
//	}
}

// returns the raw channel number for a x,y position in lab system
UInt_t TTrack::getRawChannelNumber(UInt_t det, Float_t xPred, Float_t yPred)
{
	UInt_t plane = det / 2;
	int rawChannel = 9999;
	// x planes:
	if (det%2 == 0) {
		//TODO
		rawChannel = (int)((xPred-this->getXOffset(plane) - (yPred/*-this->getYOffset(plane)*/)*TMath::Tan(this->getPhiXOffset(plane))) * TMath::Cos(this->getPhiXOffset(plane)));
	}
	// y planes:
	else {
		this->getYOffset(plane);
		this->getPhiYOffset(plane);
	}
	// telescope
	if (det < 8) {
		if (rawChannel > -1 && rawChannel < 256) {
			return rawChannel;
		}
		else return 9999;
	}
	// diamond
	else {
		if (rawChannel > -1 && rawChannel < 128) {
			return rawChannel;
		}
		else return 9999;
	}
}

Float_t TTrack::getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred)
{
	UInt_t plane = det/2;
	// x plane
	if (det%2 == 0) {
		return this->getXPositionInDetSystem(plane, xPred, yPred);
	}
	// y plane
	else {
		return this->getYPositionInDetSystem(plane, xPred, yPred);
	}
}

Float_t TTrack::getXPositionInDetSystem(UInt_t plane, Float_t xPred, Float_t yPred)
{
	return (xPred-this->getXOffset(plane))/TMath::Cos(this->getPhiXOffset(plane)) - (yPred+(xPred-this->getXOffset(plane))*TMath::Tan(this->getPhiXOffset(plane)))*TMath::Sin(this->getPhiXOffset(plane));
}

Float_t TTrack::getYPositionInDetSystem(UInt_t plane, Float_t xPred, Float_t yPred)
{
	return (yPred-this->getYOffset(plane))/TMath::Cos(this->getPhiYOffset(plane)) + (xPred-(yPred-this->getYOffset(plane))*TMath::Tan(this->getPhiYOffset(plane)))*TMath::Sin(this->getPhiYOffset(plane));
}

UInt_t TTrack::getVerbosity() const
{
    return verbosity;
}

void TTrack::setEtaIntergral(UInt_t det, TH1F *histo)
{
	cout<<"TTRack set Eta Integral of histoMap of detector "<<det<<endl;
	histoMap[det]=(TH1F*)histo->Clone();
}

TH1F *TTrack::getEtaIntergral(UInt_t det)
{
//	cout<<"get etaIntegral of det "<<det<<" in " <<histoMap.size()<<" histos:"<<flush;
	if(histoMap.find(det)!=histoMap.end() ){
//		cout<<" found histo "<<histoMap[det]->GetTitle()<<endl;
		return histoMap[det];
	}
	else
		return 0;
}








