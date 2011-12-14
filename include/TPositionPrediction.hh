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
class TPositionPrediction:public TObject {
public:
	TPositionPrediction(){xPos=0;yPos=0;xSigma=0;ySigma=0;xChi2=0;yChi2=0;bValid=false;};
	TPositionPrediction(Float_t xPos,Float_t xSigma,Float_t xChi2, Float_t yPos, Float_t ySigma, Float_t yChi2);
	virtual ~TPositionPrediction();
	Float_t getPosition(TPlane::enumCoordinate cor);
	Float_t getSigma(TPlane::enumCoordinate cor);
	Float_t getChi2(TPlane::enumCoordinate cor);
    void setxPos(Float_t pos);
    void setxSigma(Float_t sigma);
    void setxChi2(Float_t chi2);
    void setyPos(Float_t pos);
    void setySigma(Float_t sigma);
    void setyChi2(Float_t chi2);
    Float_t getPositionX(){return xPos;};
	Float_t getSigmaX(){return xSigma;};
	Float_t getChi2X(){return xChi2;};
	Float_t getPositionY(){return yPos;};
	Float_t getSigmaY(){return ySigma;};
	Float_t getChi2Y(){return yChi2;};
	bool isValid(){return bValid;};
	void setValid(bool bValid){this->bValid=bValid;};
	void Print();
private:
	Float_t xPos;
	Float_t yPos;
	Float_t xSigma;
	Float_t ySigma;
	Float_t xChi2;
	Float_t yChi2;
	bool bValid;
};

#endif /* TPOSITIONPREDICTION_HH_ */
