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
#include <vector>
#define N_INVALID -9999
class TPlaneProperties:public TObject {
public:
	enum enumCoordinate{ X_COR =0, Y_COR=1, Z_COR =2, XY_COR=3,};
	enum enumDetectorType{kUndefined = 0, kSilicon = 1, kDiamond =2};
	TPlaneProperties();
	virtual ~TPlaneProperties();
	inline static Double_t getStripDistance(){return 50.;};//strip distance in mum
	inline static UInt_t getNChannelsSilicon(){return 256;};
	inline static UInt_t getNChannelsDiamond(){return 128;};
	inline static UInt_t getNChannels(UInt_t det){switch (det){case 8: return TPlaneProperties::getNChannelsDiamond();break;default: return TPlaneProperties::getNChannelsSilicon();break;}}
	inline static UInt_t getMaxSignalHeightSilicon(){return 255;};
	inline static UInt_t getMaxSignalHeightDiamond(){return 4095;};
	inline static Int_t getMaxSignalHeight(UInt_t det){switch(det){case 8: return getMaxSignalHeightDiamond();default: return getMaxSignalHeightSilicon();}}
	inline static UInt_t getPlaneNumber(UInt_t det){return det/2;};
	inline static UInt_t getNSiliconPlanes(){return 4;};
	inline static UInt_t getNPlanes(){return getNSiliconPlanes()+1;};
	inline static UInt_t getNSiliconDetectors(){return 8;};
	inline static UInt_t getDetDiamond(){return 8;};
	inline static UInt_t getDiamondPlane(){return 4;};
	static std::vector<UInt_t> getSiliconPlaneVector();
	inline static UInt_t getNDetectors(){return 9;};
	inline static UInt_t getMaxTransparentClusterSize(UInt_t det){return 10;};
	static std::string getCoordinateString(enumCoordinate cor);
	static std::string getDetectortypeString(enumDetectorType type);
	static std::string getDetectorNameString(UInt_t det);
	static std::string getStringForDetector(int i);
	inline static bool isSiliconDetector(UInt_t det){return (det<getNSiliconDetectors());}
	inline static bool isDiamondDetector(UInt_t det){return (det==getDetDiamond());}
	inline static bool isDiamondPlane(UInt_t plane){return (plane==getDiamondPlane());}
	inline static bool isSaturated(UInt_t det, Int_t adcValue){return (adcValue>=getMaxSignalHeight(det));}

    ClassDef(TPlaneProperties,1);

};

#endif /* TPLANEPROPERTIES_HH_ */
