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
TTrack::TTrack(TDetectorAlignment *alignment,TSettings* settings) {
	verbosity=settings->getVerbosity();
	this->alignment = alignment;
	this->settings = settings;
	event=NULL;

	linFitX = new TLinearFitter(1,"pol1","D");
	linFitY= new TLinearFitter(1,"pol1","D");
	linFitX->StoreData(true);
	linFitY->StoreData(true);

	linFitX->SetFormula("pol1");
	linFitY->SetFormula("pol1");
}

TTrack::~TTrack() {
	for(std::map<UInt_t,TH1F*>::iterator iterator = histoMap.begin(); iterator != histoMap.end(); iterator++) {
//		UInt_t det = iterator->first;
		if (iterator->second) delete iterator->second;
		histoMap.erase(iterator);
	}
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

Float_t TTrack::getXMeasuredClusterPositionMetricSpace(UInt_t plane,bool cmnCorrected, TCluster::calculationMode_t mode, TH1F* histo) {
    if (verbosity>6&&TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
        cout<<"[TTrack::getXMeasuredClusterPositionMetricSpace] "<<mode<<" "<<histo<<endl;
	return this->getMeasuredClusterPositionMetricSpace(plane*2,cmnCorrected,mode,histo);
}

Float_t TTrack::getYMeasuredClusterPositionMetricSpace(UInt_t plane,bool cmnCorrected, TCluster::calculationMode_t mode, TH1F* histo) {
    if (verbosity>6&&TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
            cout<<"[TTrack::getYMeasuredClusterPositionMetricSpace] "<<mode<<" "<<histo<<endl;
    return this->getMeasuredClusterPositionMetricSpace(plane*2+1,cmnCorrected,mode,histo);
}

Float_t TTrack::getMeasuredClusterPositionMetricSpace(TPlaneProperties::enumCoordinate cor, UInt_t plane,bool cmnCorrected, TCluster::calculationMode_t mode, TH1F* histo) {
    if (verbosity>6&&TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
            cout<<"[TTrack::getMeasuredClusterPositionMetricSpace_COR_PLANE] cor"<<cor<<" plane "<<plane<<" "<<mode<<" "<<histo<<endl;
    if (cor==TPlaneProperties::X_COR)
		return getMeasuredClusterPositionMetricSpace(plane*2,cmnCorrected,mode,histo);
	else if (cor==TPlaneProperties::Y_COR)
		return getMeasuredClusterPositionMetricSpace(plane*2+1,cmnCorrected,mode,histo);
	return N_INVALID;
}

Float_t TTrack::getMeasuredClusterPositionMetricSpace(UInt_t det,bool cmnCorrected, TCluster::calculationMode_t mode, TH1F* histo) {

    if (verbosity>6&&TPlaneProperties::isSiliconDetector(det) && (mode != TCluster::corEta||histo==0))
        cout<<"[TTrack::getMeasuredClusterPositionMetricSpace_DET] det "<<det<<" "<<mode<<" "<<histo<<endl;

	return inMetricDetectorSpace(det,this->getMeasuredClusterPositionChannelSpace(det,cmnCorrected,mode,histo));
}

Float_t TTrack::inMetricDetectorSpace(UInt_t det, Float_t clusterPosition){
	Float_t metricValue = settings->convertChannelToMetric(det,clusterPosition);
	return metricValue;
}


Float_t TTrack::inChannelDetectorSpace(UInt_t det, Float_t metricPosition){
	Float_t channelPosition  = N_INVALID;
	if(TPlaneProperties::isSiliconDetector(det))
		channelPosition =  metricPosition/settings->getSiliconPitchWidth();
	else if(TPlaneProperties::isDiamondDetector(det)){
//		cout<<"validPAttern  TTrack::inChannelDetectorSpace "<<det<<" "<<metricPosition<<": "<<flush;
//		cout<<settings->diamondPattern.hasInvalidIntervals()<<endl;
		channelPosition = settings->diamondPattern.convertMetricToChannel(metricPosition);
	}
	return channelPosition;
}
/**
 * calculate Position of hit in the lab frame using the  X/Y detector cluster positions in channel numbers.
 * The way to calculate the cluster position can be selected by the mode varibale.
 *
 * @param cor coordinate in which the calculation is made
 * @param plane which plane no is used
 * @param xCluster
 * @param yCluster
 * @param mode
 * @return value of calculated position
 */
Float_t TTrack::getPostionInLabFrame(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster xCluster,TCluster yCluster,bool cmnCorrected, TCluster::calculationMode_t mode,TH1F* histo){
	if(xCluster.size()<=0||yCluster.size()<=0)
		return N_INVALID;
	if(TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
	    cerr<<"Silicon Plane "<<plane<<" but wrong mode or histo==0: "<<mode<<" "<<histo<<endl;
	Float_t xOffMeas = this->getXOffset(plane);
	Float_t yOffMeas = this->getYOffset(plane);
	Float_t xPhiOff = this->getPhiXOffset(plane);
	Float_t yPhiOff = this->getPhiYOffset(plane);
	UInt_t detX = plane*2;
	UInt_t detY = plane*2+1;
	Float_t xMeasMetric = inMetricDetectorSpace(detX,xCluster.getPosition(cmnCorrected,mode,getEtaIntegral(detX)));
	Float_t yMeasMetric = inMetricDetectorSpace(detY,yCluster.getPosition(cmnCorrected,mode,getEtaIntegral(detY)));

	// get gradients of straight lines
	Float_t m_x = xPhiOff==0?0:1./TMath::Tan(xPhiOff);//cotan
	Float_t m_y = TMath::Tan(-1.*yPhiOff);//tan
	//take Position in the detecor frame at y_det = 0 and convert it into the lab frame
	Float_t x_x = xOffMeas + TMath::Cos(xPhiOff)*xMeasMetric;
	Float_t y_x = (-1.)*TMath::Sin(xPhiOff)*xMeasMetric;
	//take y Position in the dector frame at x_det = 0 and convert it into the lab frame
	Float_t x_y = TMath::Sin(yPhiOff)*yMeasMetric;
	Float_t y_y = yOffMeas+TMath::Cos(yPhiOff)*yMeasMetric;

	//out of slope and position we can calculate the offset
	Float_t b_x = xPhiOff==0?xOffMeas+xMeasMetric:y_x-m_x*x_x;
	Float_t b_y = y_y-m_y*x_y;
	//find interception of both straight liness
	Float_t x = (b_y-b_x)/(m_x-m_y);
	if (xPhiOff==0) x = xOffMeas+xMeasMetric;
	Float_t xPosition = x;
	Float_t yPosition = m_y*x+b_y;
	if(m_y!=m_y){
		cerr<<"There was a problem with m_y! Please have a look. yPhiOff: "<<yPhiOff<<endl;
		exit(-1);
	}
	if(m_x!=m_x){
		cerr<<"There was a problem with m_y! Please have a look. xPhiOff: "<<xPhiOff<<endl;
		exit(-1);
	}

	if(verbosity&&(xPosition==-1||yPosition==-1)){
		cout<<" Meas: "<<xMeasMetric<<"/"<<yMeasMetric<<"\toff: "<<xOffMeas<<"/"<<yOffMeas<<"\tphi"<<xPhiOff<<"/"<<yPhiOff<<endl;
		cout<<"X:"<<mode<<" "<<xCluster.getPosition(cmnCorrected,mode,histo)<<endl;
		xCluster.Print(2);
		cout<<"Y: "<<endl;
		yCluster.Print(2);
	}

	switch(cor){
		case TPlaneProperties::X_COR: return xPosition;break;
		case TPlaneProperties::Y_COR: return yPosition;break;
		case TPlaneProperties::Z_COR: return getZPosition(plane);break;
		default: return N_INVALID;
	}
}

/**
 * calculate Position in lab system of a cluster when there is no perpendicular strip detector
 *
 * @param det detector in which the cluster is
 * @param cluster
 * @param predictedPerpPosition predicted position of the track in the perpendicular det of the same plane (in lab system)
 * @param mode
 * @return value of calculated position
 */
Float_t TTrack::getPositionOfCluster(UInt_t det, TCluster cluster, Float_t predictedPerpPosition, bool cmnCorrected, TCluster::calculationMode_t mode,TH1F* histo) {
	if(mode == TCluster::corEta && histo==0)
		cout<<"getPositionOfCluster::histo ==0"<<endl;
	if (cluster.size()<=0) return N_INVALID;
	UInt_t plane = det/2;
	if(histo==0)
		histo = (TH1F*)getEtaIntegral(plane*2+1);
	Float_t meas = inMetricDetectorSpace(det,cluster.getPosition(cmnCorrected,mode,histo));
	Float_t xOffMeas = this->getXOffset(plane);
	Float_t yOffMeas = this->getYOffset(plane);
	Float_t xPhiOff = this->getPhiXOffset(plane);
	Float_t yPhiOff = this->getPhiYOffset(plane);
	Float_t m_x = 1./TMath::Tan(xPhiOff);
	Float_t m_y = (-1.)*TMath::Tan(yPhiOff);
	Float_t x_x = xOffMeas + TMath::Cos(xPhiOff)*meas;
	Float_t y_x = (-1.)*TMath::Sin(xPhiOff)*meas;
	Float_t x_y = TMath::Sin(yPhiOff)*meas;
	Float_t y_y = yOffMeas+TMath::Cos(yPhiOff)*meas;
	Float_t b_x = y_x-m_x*x_x;
	Float_t b_y = y_y-m_y*x_y;
	Float_t xPosition = (predictedPerpPosition-b_x)/m_x;
	Float_t yPosition = m_y*predictedPerpPosition+b_y;
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
 * @return calculated hit position in Metric Space
 */
Float_t TTrack::getStripXPosition(UInt_t plane,Float_t yPred,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo){
	if(event==NULL)return N_INVALID;
	if(event->getNXClusters(plane)!=1)
		return N_INVALID;
	// get offsets
	TCluster xCluster = event->getCluster(plane,TPlaneProperties::X_COR,0);
	if(xCluster.size()<=0) {
		cerr<<"This cluster is too small!!!!"<<endl;
		xCluster.Print();
	}
	if(histo==NULL)
		histo = getEtaIntegral(plane*2);
	return getXPositionInLabFrameStripDetector(plane,xCluster,yPred,cmnCorrected, mode,histo);
}

/**
 * Calculation of Hitposition of event using the measured Offsets
 * in this transformation it transform the measured hitPosition in the
 * plane space into the space of the first plane
 * @param cor enum of Coordinate from which you want to have the result
 * @param plane number of plane for which you are calculating the hitposition
 * @param yPred prediction of Y Position in metric Space
 * @return calculated hit position in Metric LAB Frame
 */
Float_t TTrack::getXPositionInLabFrameStripDetector(UInt_t plane,TCluster xCluster, Float_t yPred,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo){
	if(xCluster.size()<=0)
		return N_INVALID;
	// get offsets
	UInt_t det = plane * 2;
	Float_t xOffset = this->getXOffset(plane);
	Float_t phiXOffset = this->getPhiXOffset(plane);
	if(histo==NULL)
		histo = getEtaIntegral(det);
	Float_t xMeasured = inMetricDetectorSpace(det,xCluster.getPosition(cmnCorrected,mode,histo));//-xOffset;
	// apply offsets
	//calculation checked. is correct
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
Float_t TTrack::getXPositionMetric(UInt_t plane,bool cmnCorrected, TCluster::calculationMode_t mode,TH1F* histo) {
	if(histo==0)
		histo = getEtaIntegral(plane*2);
	return getPositionInLabFrame(TPlaneProperties::X_COR,plane,cmnCorrected,mode,histo);
}

/**
 * gives the yPosition of "event" in the telescope space
 * this function uses getPosition(Y_COR,plane)
 * @param plane planeNumber to use correct offsets
 * @return calculated yPosition
 */
Float_t TTrack::getYPositionMetric(UInt_t plane,bool cmnCorrected, TCluster::calculationMode_t mode,TH1F* histo) {
	if(histo==0)
			histo = getEtaIntegral(plane*2);
	return getPositionInLabFrame(TPlaneProperties::Y_COR,plane,cmnCorrected,mode,histo);
}

Float_t TTrack::getZPosition(UInt_t plane,TCluster::calculationMode_t mode){
	return alignment->GetZOffset(plane);
}


Float_t TTrack::calculateTrackResolution(Float_t zPosSubjectPlane,
		vector<Float_t> zPositions, vector<Float_t> resolutions) {
//	cout<<"TTrack::calculateTrackResolution "<<zPositions.size()<<" "<<resolutions.size()<<endl;
//	Float_t zSigma= alignment->getZResolution(subjectPlane);//todo
	linFitX->ClearPoints();
	for(UInt_t pl=0;pl<zPositions.size();pl++){
		vector<Double_t> zPosVec;
//		cout<<"adding "<<pl<<"\t"<<flush;
		Float_t x = zPositions.at(pl)-zPosSubjectPlane;
		Float_t y = 0;
		Float_t res = resolutions.at(pl);
//		cout<<"adding "<<pl<<" "<<x<<" "<<y<<" "<<res<<flush;
		zPosVec.push_back(x);
//		cout<<"#"<<endl;
		linFitX->AddPoint(&zPosVec.at(0),y,res);
//		cout<<"DONE"<<endl;
	}
//	cout<<"evaluate"<<endl;
	linFitX->Eval();
//	cout<<"evaluated"<<linFitX<<"\nGetParameter 1"<<endl;
	Float_t m = linFitX->GetParameter(1);
//	cout <<"Get ParError 1"<<endl;
	Float_t sigma_m = linFitX->GetParError(1);
//	cout<<"GetParameter 0"<<endl;
	Float_t b = linFitX->GetParameter(0);
//	cout<<"Get PAr Error 0"<<endl;
	Float_t sigma_b = linFitX->GetParError(0);
//	cout<<"get result"<<endl;
	pair<Float_t,Float_t> result = getPredictedHitPosition(0,m,b,0,sigma_m,sigma_b);
	return result.second;
}


std::pair<Float_t, Float_t> TTrack::getPredictedHitPosition(Float_t zPos,Float_t m, Float_t b,Float_t sigma_z, Float_t sigma_m, Float_t sigma_b){
	Float_t position = m*zPos+b;
	Float_t sigma = (zPos*sigma_m)*(zPos*sigma_m)+(m*sigma_z)*(m*sigma_z)+(sigma_b*sigma_b);
	sigma = TMath::Sqrt(sigma);
	return make_pair(position,sigma);
}
/**
 * @brief predicts the Position calculated by a linear fit out of all planes in vecRefPlanes
 * In Order to be able to calculate the right sigma the subject Plane is set to z_subject = 0
 * all other distances to the planes are set to the difference.
 * This must be done in order to calculate the sigma correctly:
 * sigma_b and sigma_m are correlated to each other in order to make it easy as possible we set the
 * subjectplane to z=0 such that sigma_prediction**2 = sigma_b**2 since the term sigma_m**2*zPos**2
 * is equal to zero
 * @param subjectPlane Plane where the Position should be predicted
 * @param vecRefPlanes Planes used for the prediction
 * @param mode the way how the position should be calculated, e.g. Eta corrected
 * @param bPrint show output?
 * @return PositionPrediction object which contains all Informations of the prediction.
 */
TPositionPrediction* TTrack::predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool cmnCorrection, TCluster::calculationMode_t mode,bool bPrint)
{
    if(mode!=TCluster::corEta)
            cout<<"[TTrack::predictPosition] strange mode: "<<mode<<endl;
	Float_t zPosSubjectPlane = getZPosition(subjectPlane);
	Float_t sigma_z= alignment->getZResolution(subjectPlane);//todo
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
		if(verbosity>10||bPrint)	cout<<"TTrack::predictPosition with 1 refPlane"<<endl;
		//todo anpassen so dass sigmas da drin reinkommen...
		TPositionPrediction *prediction=new TPositionPrediction(zPosSubjectPlane,sigma_z,
		        getXPositionMetric(vecRefPlanes.at(0),cmnCorrection,mode), 0.,0.,
		        getYPositionMetric(vecRefPlanes.at(0),cmnCorrection,mode),0.,0.,0,0);
		return prediction;
	}
	vector<Double_t> zPosVec;//todo add xsigma ysigma
	if(bPrint)cout<<"Prediction of Track in Plane "<<subjectPlane<<" with "<<vecRefPlanes.size()<<" Planes:"<<endl;
	bool lastPredictionValid=true;
	for(UInt_t pl=0;pl<vecRefPlanes.size();pl++){
		UInt_t plane=vecRefPlanes.at(pl);
		zPosVec.clear();
		Float_t zPos = alignment->GetZOffset(plane)-zPosSubjectPlane;
		zPosVec.push_back(zPos);
		Float_t xPos = (Double_t)getXPositionMetric(plane,cmnCorrection,mode,getEtaIntegral(pl*2));
		Float_t yPos = (Double_t)getYPositionMetric(plane,cmnCorrection,mode,getEtaIntegral(pl*2+1));
		if((xPos==-1||yPos==-1)&&(verbosity||bPrint)){
		    if(bPrint)cout<<"Problem with Plane "<<plane<<" "<<xPos<<" "<<yPos<<endl;
			event->Print(1);
		}
		Float_t xRes = this->alignment->getXResolution(plane);
		if (xRes <= 0)  {
		    xRes = 10;
		    if(bPrint)cout<<"UPdate xREs to "<<xRes<<endl;
		}
		Float_t yRes = this->alignment->getYResolution(plane);
		if (yRes <= 0)  {
		    yRes = 10;
		    if(bPrint)cout<<"Update yRes to "<<yRes<<endl;
		}
		linFitX->AddPoint(&zPosVec.at(0),xPos,xRes);//todo anpassen des SIGMA
		linFitY->AddPoint(&zPosVec.at(0),yPos,yRes);//todo anpassen des sigma 0.001
		if(xPos==-1||yPos==-1)
			lastPredictionValid = false;
		if(verbosity>10||bPrint)
			cout<<"\tAdd in Plane "<<plane<<"  X: "<<xPos<<"+/-"<<xRes<<" / "<<
			                                "  Y: "<<yPos<<"+/-"<<yRes<<" / "<<
			                                "  Z: "<<getZPosition(plane) <<"("<<zPosVec.at(0)<<")"<<endl;
	}
	linFitX->Eval();
	linFitY->Eval();
	linFitX->Chisquare();
	linFitY->Chisquare();
//	pair<Float_t,Float_t> result = getPredictedHitPosition(zPos,mx,bx,sigma_z,sigma_mx,sigma_bx)
	Float_t zPos = 0;//zPosSubjectPlane;//alignment->GetZOffset(subjectPlane);
	Float_t mx = linFitX->GetParameter(1);
	Float_t sigma_mx = linFitX->GetParError(1);
	Float_t bx = linFitX->GetParameter(0);
	Float_t sigma_bx = linFitX->GetParError(0);
	Float_t my = linFitY->GetParameter(1);
	Float_t sigma_my = linFitY->GetParError(1);
	Float_t by = linFitY->GetParameter(0);
	Float_t sigma_by = linFitY->GetParError(0);
	Float_t xChi2 = linFitX->GetChisquare()/(linFitX->GetNpoints()-linFitX->GetNumberFreeParameters());
	Float_t yChi2 = linFitY->GetChisquare()/(linFitY->GetNpoints()-linFitY->GetNumberFreeParameters());
	if(verbosity>10){
			cout<<"\tParameters:\n\t mx: "<<mx<<" +/- "<<sigma_mx<<"\tbx: "<<bx<<" +/- "<<sigma_bx<<"\tzPos:"<<zPos<<endl;
			cout<<"\t my: "<<my<<" +/- "<<sigma_my<<"\tby: "<<by<<" +/- "<<sigma_by<<"\tzPos:"<<zPos<<endl;
			cout<<"\tNDFx: "<<linFitX->GetNumberFreeParameters()<<"\t"<<"\tNDFy: "<<linFitY->GetNumberFreeParameters()<<endl;
			cout<<"\t chi2x:"<<xChi2<<"\tchi2y:"<<yChi2<<endl;
		}
	pair<Float_t,Float_t> xprediction = this->getPredictedHitPosition(zPos,mx,bx,sigma_z,sigma_mx,sigma_bx);
	pair<Float_t,Float_t> yprediction = this->getPredictedHitPosition(zPos,my,by,sigma_z,sigma_my,sigma_by);
//	Float_t xPos = mx*zPos+bx;
//	Float_t yPos = my*zPos+by;
//	Float_t xSigma = (zPos*sigma_mx)*(zPos*sigma_mx)+(mx*sigma_z)*(mx*sigma_z)+(sigma_bx*sigma_bx);
//	xSigma = TMath::Sqrt(xSigma);
//	Float_t ySigma = (zPos*sigma_my)*(zPos*sigma_my)+(my*sigma_z)*(my*sigma_z)+(sigma_by*sigma_by);
//	ySigma = TMath::Sqrt(ySigma);
	Float_t xPhi = TMath::ATan(mx);
	Float_t yPhi = TMath::ATan(my);
	zPos=zPosSubjectPlane;
	TPositionPrediction* prediction=new TPositionPrediction(zPos,sigma_z,xprediction.first,xprediction.second,xChi2,yprediction.first,yprediction.second,yChi2,xPhi,yPhi);
	if(!lastPredictionValid)prediction->setValid(false);
	if(verbosity>10||bPrint)	cout<<"\tPrediction of Plane "<<subjectPlane<<" with "<<vecRefPlanes.size()<<" Planes for ZPosition: "<<zPos<<endl;
	if(verbosity>10||bPrint)	cout<<"\t X: "<<xprediction.first<<" +/- "<<xprediction.second<<"   with a Chi^2 of "<<xChi2<<"  "<<linFitX->GetNpoints()<<endl;
	if(verbosity>10||bPrint)	cout<<"\t Y: "<<yprediction.first<<" +/- "<<yprediction.second<<"   with a Chi^2 of "<<yChi2<<"  "<<linFitY->GetNpoints()<<"\n";
    if(verbosity>10||bPrint)	cout<<"\t PhiX: "<<xPhi<<"  PhiY:"<<yPhi<<"\tm:"<<mx<<"/"<<my<<"  b:"<<bx<<"/"<<by<<"\n"<<endl;
	return prediction;
}

vector<Float_t> TTrack::getSiXPositions() {
	vector<Float_t> xPositions;
	if(event==NULL)return xPositions;
	for (int plane = 0; plane < 4; plane++) {
		xPositions.push_back(this->getXPositionMetric(plane,false));
	}
	return xPositions;
}

vector<Float_t> TTrack::getSiYPositions() {

	vector<Float_t> yPositions;
	if(event==NULL)return yPositions;
	for (int plane=0; plane < 4; plane++) {
		yPositions.push_back(this->getYPositionMetric(plane,false));
	}
	return yPositions;
}

void TTrack::setDetectorAlignment(TDetectorAlignment *alignment)
{
	if(alignment!=this->alignment)
		delete this->alignment;
	this->alignment=alignment;
}

Float_t TTrack::getXMeasuredClusterPositionChannelSpace(UInt_t plane,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo)
{
    if (verbosity>6&&TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
            cout<<"[TTrack::getXMeasuredClusterPositionChannelSpace] "<<mode<<" "<<histo<<endl;

	return getMeasuredClusterPositionChannelSpace
			(TPlaneProperties::X_COR,plane,cmnCorrected,mode,histo);
}

Float_t TTrack::getYMeasuredClusterPositionChannelSpace(UInt_t plane,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo)
{
    if (verbosity>6&&TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
            cout<<"[TTrack::getYMeasuredClusterPositionChannelSpace] "<<mode<<" "<<histo<<endl;
	return getMeasuredClusterPositionChannelSpace(TPlaneProperties::Y_COR,plane,cmnCorrected,mode,histo);
}
/**
 * Calculation of Hitposition of event using the measured Offsets
 * in this transformation it transform the measured hitPosition in the
 * plane space into the space of the first plane
 * @param cor enum of Coordinate from which you want to have the result
 * @param plane number of plane for which you are calculating the hitposition
 * @return calculated hit position
 */
Float_t TTrack::getPositionInLabFrame(TPlaneProperties::enumCoordinate cor,UInt_t plane, bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo){
	if(event==NULL)return N_INVALID;
	if(event->getNXClusters(plane)!=1||event->getNYClusters(plane)!=1){
		if(verbosity>5) cout<<event->getEventNumber()<<" TTRack::getPositionInLabFrame: "<<event->getNXClusters(plane)<<" "<<event->getNYClusters(plane)<<" "<<endl;
		return N_INVALID;
	}
	// get offsets
	TCluster xCluster,yCluster;
	xCluster=event->getPlane(plane).getXCluster(0);
	yCluster=event->getPlane(plane).getYCluster(0);
	if(verbosity>5)
		cout<<event<<":"<<event->getEventNumber()<<"->"<<xCluster.getPosition(cmnCorrected)<<"/"<<yCluster.getPosition(cmnCorrected)<<endl;
	return getPostionInLabFrame(cor,plane,xCluster,yCluster, cmnCorrected,mode,histo);
}


/**
 * if there exist one and only one Cluster in the plane for a given coordinate it returns the
 * channel position of the cluster calculate with a certain mode
 * @param cor
 * @param plane
 * @param mode
 * @param histo
 * @return
 */
Float_t TTrack::getMeasuredClusterPositionChannelSpace(TPlaneProperties::enumCoordinate cor, UInt_t plane,bool cmnCorrected, TCluster::calculationMode_t mode,TH1F* histo)
{
	if(event==NULL)
		return N_INVALID;
	if(cor==TPlaneProperties::XY_COR&&(event->getNXClusters(plane)!=1||event->getNYClusters(plane)!=1))
		return N_INVALID;
	if (verbosity>6&&TPlaneProperties::isSiliconPlane(plane) && (mode != TCluster::corEta||histo==0))
	    cout<<"[TTrack::getMeasuredClusterPositionChannelSpace] COR/Plane "<<mode<<" "<<histo<<endl;
	if(cor==TPlaneProperties::X_COR&&event->getNXClusters(plane)==1)
		return event->getPlane(plane).getXPosition(0,cmnCorrected,mode,histo);
	if(cor==TPlaneProperties::Y_COR&&event->getNYClusters(plane)==1)
		return event->getPlane(plane).getYPosition(0,cmnCorrected,mode,histo);
	return N_INVALID;
//// get offsets
//	switch(cor){
//	case TPlaneProperties::X_COR:break;
//	case TPlaneProperties::Y_COR:break;
//	default: return N_INVALID;
//	}
}

Float_t TTrack::getMeasuredClusterPositionChannelSpace(UInt_t det,bool cmnCorrected,TCluster::calculationMode_t mode,TH1F* histo)
{
    if (verbosity>6&&TPlaneProperties::isSiliconDetector(det) && (mode != TCluster::corEta||histo==0))
        cout<<"[TTrack::getMeasuredClusterPositionChannelSpace] DET"<<det<<" "<<mode<<" "<<histo<<endl;
	if (det%2==0)
		return getMeasuredClusterPositionChannelSpace(TPlaneProperties::X_COR,det/2,cmnCorrected,mode,histo);
	else
		return getMeasuredClusterPositionChannelSpace(TPlaneProperties::Y_COR,det/2,cmnCorrected,mode,histo);
	return N_INVALID;
}

/**
 * from a given prediction xPred/yPred this function calculates the corresponding raw channel no in det Detector
 * where the hit should have been
 * @author Lukas Baeni
 * @param det
 * @param xPred
 * @param yPred
 * @return
 * @todo Felix: I just calculated the things by my own but I'm not sure if there is an error:
 * 			with
 * 			 x_m: measured x pos (Detecotr frame)
 * 			 y_i: detector frame y position (not measured)
 * 			 x_off: offset x
 * 			 phi: angle of lab frame to detector frame
 * 			x_P/y_P = (x_off + cos phi * x_m + sin phi * y_i/ -sin phi x_m + cos phi * y_i)
 * 			solve it to get x_m:
 * 			x_m = cos phi * x(x
 */
UInt_t TTrack::getRawChannelNumber(UInt_t det, Float_t xPred, Float_t yPred)
{
	UInt_t plane = det / 2;
	int rawChannel = 9999;
	Float_t metricDetectorPosition = N_INVALID;
	// x planes:
	if (det%2 == 0) {
		//TODO
		metricDetectorPosition = (int)((xPred-this->getXOffset(plane) - (yPred/*-this->getYOffset(plane)*/)*TMath::Tan(this->getPhiXOffset(plane))) * TMath::Cos(this->getPhiXOffset(plane)));
	}
	// y planes:
	/// @todo einfuegen
	else {
		this->getYOffset(plane);
		this->getPhiYOffset(plane);
	}
	rawChannel = inChannelDetectorSpace(det,metricDetectorPosition);
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

/** from pediction calculated hit Position of a detector in detector system, metric space
 *
 * @param det
 * @param xPred
 * @param yPred
 * @return
 */
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

/** from prediction calculated X Position in Detector System - Metric Space
 *
 * @param plane
 * @param xPred
 * @param yPred
 * @return
 */
Float_t TTrack::getXPositionInDetSystem(UInt_t plane, Float_t xPred, Float_t yPred)
{
	return (xPred-this->getXOffset(plane))/TMath::Cos(this->getPhiXOffset(plane)) - (yPred+(xPred-this->getXOffset(plane))*TMath::Tan(this->getPhiXOffset(plane)))*TMath::Sin(this->getPhiXOffset(plane));
}

/** from prediction calculated Y Position in Detector System - Metric Space
 *
 * @param plane
 * @param xPred
 * @param yPred
 * @return
 */
Float_t TTrack::getYPositionInDetSystem(UInt_t plane, Float_t xPred, Float_t yPred)
{
	return (yPred-this->getYOffset(plane))/TMath::Cos(this->getPhiYOffset(plane)) + (xPred-(yPred-this->getYOffset(plane))*TMath::Tan(this->getPhiYOffset(plane)))*TMath::Sin(this->getPhiYOffset(plane));
}


/** from prediction calulation X Position in Strip Detector
 * @todo make it work for a y strip detector as well
 * @param plane
 * @param xPred
 * @param yPred
 * @return
 */
Float_t  TTrack::getXPositionInStripDetSystem(UInt_t plane, Float_t xPred, Float_t yPred){

}

/** from prediction calulation Y Position in Strip Detector
 *
 * @param plane
 * @param xPred
 * @param yPred
 * @return
 */
Float_t  TTrack::getYPositionInStripDetSystem(UInt_t plane, Float_t xPred, Float_t yPred){

	Float_t xPredOff = xPred - this->getXOffset(plane);
	Float_t phiX = getPhiXOffset(plane);
	return yPred * TMath::Cos(phiX) + xPredOff * TMath::Sin(phiX);
}


UInt_t TTrack::getVerbosity() const
{
    return verbosity;
}

void TTrack::setEtaIntegral(UInt_t det, TH1F *histo)
{

//	cout<<"TTRack set Eta Integral of histoMap of detector "<<det<<endl;
	if(histo!=0)
	  histoMap[det]=(TH1F*)histo->Clone();
	else
	  histoMap[det]=0;
}

TH1F *TTrack::getEtaIntegral(UInt_t det)
{
//	cout<<"get etaIntegral of det "<<det<<" in " <<histoMap.size()<<" histos:"<<flush;
	if(histoMap.find(det)!=histoMap.end() ){
//		cout<<" found histo "<<histoMap[det]->GetTitle()<<endl;
		return histoMap[det];
	}
	else
		return 0;
}

Float_t TTrack::getRelativeHitPosition(UInt_t det, Float_t hitPosDetMetric) {
    if(TPlaneProperties::isSiliconDetector(det))
        return fmod(hitPosDetMetric,settings->getSiliconPitchWidth());
    else if(TPlaneProperties::isDiamondDetector(det)){
//      cout<<"validPAttern  TTrack::inChannelDetectorSpace "<<det<<" "<<metricPosition<<": "<<flush;
//      cout<<settings->diamondPattern.hasInvalidIntervals()<<endl;
      Float_t relMetric = settings->diamondPattern.convertMetricToRelativeMetric(hitPosDetMetric);
      return relMetric;
    }
    return N_INVALID;
}
