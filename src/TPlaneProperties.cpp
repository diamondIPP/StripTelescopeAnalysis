/*
 * TPlaneProperties.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: bachmair
 */

#include "../include/TPlaneProperties.hh"

TPlaneProperties::TPlaneProperties() {
	// TODO Auto-generated constructor stub

}

TPlaneProperties::~TPlaneProperties() {
	// TODO Auto-generated destructor stub
}

 UInt_t TPlaneProperties::getNChannels(UInt_t det){
	switch (det){
	case 8: return TPlaneProperties::getNChannelsDiamond();break;
	default: return TPlaneProperties::getNChannelsSilicon();break;
	}
}
 UInt_t TPlaneProperties::getMaxSignalHeight(UInt_t det){
	switch(det){
	case 8: return getMaxSignalHeightDiamond();
	default: return getMaxSignalHeightSilicon();
	}
}
