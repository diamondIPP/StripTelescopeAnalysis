/*
 * TPlaneProperties.h
 *
 *  Created on: Feb 2, 2012
 *      Author: bachmair
 */

#ifndef TPLANEPROPERTIES_HH_
#define TPLANEPROPERTIES_HH_

#include <string>
#include "TROOT.h"
#include <fstream>
#include <iostream>
#include <iostream>
#include <iomanip>
#include <sstream>
class TPlaneProperties:public TObject {
public:
	enum enumCoordinate{ X_COR =0, Y_COR=1, Z_COR =2, XY_COR=3,};
	enum enumDetectorType{kUndefined = 0, kSilicon = 1, kDiamond =2};
	TPlaneProperties();
	virtual ~TPlaneProperties();
	static Double_t getStripDistance(){return 50.;};//strip distance in mum
	static UInt_t getNChannelsSilicon(){return 256;};
	static UInt_t getNChannelsDiamond(){return 128;};
	static UInt_t getNChannels(UInt_t det);
	static UInt_t getMaxSignalHeightSilicon(){return 255;};
	static UInt_t getMaxSignalHeightDiamond(){return 4095;};
	static Int_t getMaxSignalHeight(UInt_t det);
	static UInt_t getPlaneNumber(UInt_t det){return det/2;};
	static UInt_t getNSiliconPlanes(){return 4;};
	static UInt_t getNSiliconDetectors(){return 8;};
	static UInt_t getDetDiamond(){return 8;};
	static UInt_t getDiamondPlane(){return 4;};
	static UInt_t getNDetectors(){return 9;};
	static UInt_t getMaxTransparentClusterSize(UInt_t det){return 10;};
	static std::string getCoordinateString(enumCoordinate cor);
	static std::string getDetectortypeString(enumDetectorType type);
	static std::string getDetectorNameString(UInt_t det);
	static std::string getStringForDetector(int i);
	static bool isSiliconDetector(UInt_t det){return (det<getNSiliconDetectors());}
	static bool isDiamondDetector(UInt_t det){return (det==getDetDiamond());}
	static bool isDiamondPlane(UInt_t plane){return (plane==getDiamondPlane());}
	static bool isSaturated(UInt_t det, Int_t adcValue){return (adcValue>=getMaxSignalHeight(det));}

    ClassDef(TPlaneProperties,1);

};

#endif /* TPLANEPROPERTIES_HH_ */
