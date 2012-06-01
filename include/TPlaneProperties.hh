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
	static UInt_t getNChannelsSilicon(){return 256;};
	static UInt_t getNChannelsDiamond(){return 128;};
	static UInt_t getNChannels(UInt_t det);
	static UInt_t getMaxSignalHeightSilicon(){return 255;};
	static UInt_t getMaxSignalHeightDiamond(){return 4095;};
	static UInt_t getMaxSignalHeight(UInt_t det);
	static UInt_t getPlaneNumber(UInt_t det){return det/2;};
	static UInt_t getNSiliconPlanes(){return 4;};
	static UInt_t getNSiliconDetectors(){return 8;};
	static UInt_t getDetDiamond(){return 8;};
	static UInt_t getDiamondPlane(){return 4;};
	static UInt_t getNDetectors(){return 9;};
	static UInt_t getMaxTransparentClusterSize(UInt_t det){return 10;};
    static std::string getCoordinateString(enumCoordinate cor);
    static std::string getDetectortypeString(enumDetectorType type);//todo verschieben
    static std::string getDetectorNameString(UInt_t det);
	// TODO: getPlaneCoordinate(UInt_t plane);

    ClassDef(TPlaneProperties,1);

};

#endif /* TPLANEPROPERTIES_HH_ */
