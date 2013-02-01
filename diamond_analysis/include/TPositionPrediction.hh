/*
 * TPositionPrediction.hh
 *
 *  Created on: Dec 14, 2011
 *      Author: bachmair
 */

#ifndef TPOSITIONPREDICTION_HH_
#define TPOSITIONPREDICTION_HH_

#include "TPlane.hh"
#include "TObject.h"
#include "TCluster.hh"
class TPositionPrediction:public TObject {
public:
	TPositionPrediction(){xPos=0;yPos=0;xSigma=0;ySigma=0;xChi2=0;yChi2=0;bValid=false;xPhi=0;yPhi=0;zSigma = 0;zPos=0;};
	TPositionPrediction(Float_t zPos,Float_t zSigma, Float_t xPos,Float_t xSigma,Float_t xChi2, Float_t yPos, Float_t ySigma, Float_t yChi2,Float_t xPhi,Float_t yPhi);

	virtual ~TPositionPrediction();
	Float_t getPosition(TPlaneProperties::enumCoordinate cor);
	Float_t getSigma(TPlaneProperties::enumCoordinate cor);
	Float_t getChi2();
	Float_t getChi2(TPlaneProperties::enumCoordinate cor);
	Float_t getPhi(TPlaneProperties::enumCoordinate cor);
	void setxPos(Float_t pos);
	void setxSigma(Float_t sigma);
	void setxChi2(Float_t chi2);
	void setyPos(Float_t pos);
	void setySigma(Float_t sigma);
	void setyChi2(Float_t chi2);
	/**
	 * Predicted X Position in Metric Space
	 * @return
	 */
	Float_t getPositionX(){return xPos;if(xPos>-1000&&xPos<100000)return xPos;else return -9999;};
	Float_t getSigmaX(){return xSigma;};
	Float_t getChi2X(){return xChi2;};
	/**
	 * Predicted Y Position in Metric Space
	 * @return
	 */
	Float_t getPositionY(){return yPos;if(yPos>-1000&&yPos<100000)return yPos;else return -9999;};
	Float_t getSigmaY(){return ySigma;};
	Float_t getChi2Y(){return yChi2;};
	Float_t getPositionZ(){return zPos;if(zPos>-1000&&zPos<100000)return zPos;else return -9999;};
	Float_t getSigmaZ(){return zSigma;};
	Float_t getPhiX(){return xPhi;};
	Float_t getPhiY(){return yPhi;};;
	bool isValid(){return bValid;};
	void setValid(bool bValid){this->bValid=bValid;};
	void Print(UInt_t level=0);
private:
	Float_t xPos;
	Float_t yPos;
	Float_t xSigma;
	Float_t ySigma;
	Float_t xChi2;
	Float_t yChi2;
	Float_t zPos;
	Float_t zSigma;
	Float_t xPhi,yPhi;
	bool bValid;
};

#endif /* TPOSITIONPREDICTION_HH_ */
