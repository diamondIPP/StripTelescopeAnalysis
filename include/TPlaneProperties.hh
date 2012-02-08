/*
 * TPlaneProperties.h
 *
 *  Created on: Feb 2, 2012
 *      Author: bachmair
 */

#ifndef TPLANEPROPERTIES_HH_
#define TPLANEPROPERTIES_HH_

#include "TROOT.h"
class TPlaneProperties {
public:

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
	// TODO: getPlaneCoordinate(UInt_t plane);

};

#endif /* TPLANEPROPERTIES_HH_ */
