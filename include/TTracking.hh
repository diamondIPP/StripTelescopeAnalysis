/*
 * TTracking.hh
 *
 *  Created on: Jan 23, 2012
 *      Author: bachmair
 */

#ifndef TTRACKING_HH_
#define TTRACKING_HH_

#include "TEventReader.hh"
#include "TTrack.hh"
#include "TDetectorAlignment.hh"
#include "TFile.h"
#include "TPositionPrediction.hh"
#include "TSettings.class.hh"
class TTracking: public TADCEventReader{
public:
	TTracking(std::string pathName, std::string alignmentName,std::string etaDistributionPath, TSettings* settings);
	TPositionPrediction* predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bPrint=false);
	Float_t getXPosition(UInt_t plane);
	Float_t getYPosition(UInt_t plane);
	Float_t getZPosition(UInt_t plane);
	virtual ~TTracking();
	bool LoadEvent(UInt_t eventNumber);
	Float_t getStripXPositionOfCluster(UInt_t plane,TCluster xCluster, Float_t yPred,TCluster::calculationMode_t mode=TCluster::highest2Centroid,TH1F* histo=NULL);
	Float_t getStripXPosition(UInt_t plane,Float_t yPred,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getPositionOfCluster(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster xCluster,TCluster yCluster, TCluster::calculationMode_t mode=TCluster::highest2Centroid, TH1F* histo=0);
	Float_t getPositionOfCluster(UInt_t det, TCluster cluster, Float_t predictedPerpPosition, TCluster::calculationMode_t mode=TCluster::highest2Centroid, TH1F* histo=0);
	Float_t getPosition(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
	Float_t getPositionInDetSystem(UInt_t det, Float_t xPred, Float_t yPred);
	Float_t getMeasuredPositionMetricSpace(TPlaneProperties::enumCoordinate cor, UInt_t plane,TCluster::calculationMode_t mode=TCluster::highest2Centroid);
//	Float_t getPosition(TPlaneProperties::enumCoordinate cor,UInt_t plane,TCluster::calculationMode_t mode=TCluster::highest2Centroid){return myTrack->getPosition(cor,plane, mode);};
//	Float_t getXPosition(UInt_t plane);
//	Float_t getYPosition(UInt_t plane);
//	Float_t getZPosition(UInt_t plane);
//	Float_t getXMeasured(UInt_t plane);
//	Float_t getYMeasured(UInt_t plane);
//	Float_t getMeasured(TPlaneProperties::enumCoordinate cor,UInt_t plane);
//	TPositionPrediction* predictPosition(UInt_t subjectPlane, vector<UInt_t> vecRefPlanes,bool bPrint=false);
	TDetectorAlignment* getAlignment(){return myAlignment;};
private:
	bool setAlignment(std::string alignmentName);
	TTrack *myTrack;
	TSettings* settings;
	TFile* alignmentFile;
	TDetectorAlignment* myAlignment;
};

#endif /* TTRACKING_HH_ */
