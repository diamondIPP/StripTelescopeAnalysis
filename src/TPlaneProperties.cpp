/*
 * TPlaneProperties.cpp
 *
 *  Created on: Feb 2, 2012
 *      Author: bachmair
 */

#include "../include/TPlaneProperties.hh"
ClassImp(TPlaneProperties);
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

std::string TPlaneProperties::getCoordinateString(enumCoordinate cor){
	switch (cor){
	case X_COR: return "X";break;
	case Y_COR: return "Y";break;
	case Z_COR: return "Z";break;
	case XY_COR:return "X&Y"; break;
	default: return "UNDEFINDED";
	}
}

std::string TPlaneProperties::getDetectortypeString(TPlaneProperties::enumDetectorType type){
	switch (type){
	case TPlaneProperties::kSilicon: 	return "Silicon";
	case TPlaneProperties::kDiamond:	return "Diamond";
	default:		return "UNDEFINED";
	}
}


std::string TPlaneProperties::getDetectorNameString(UInt_t det){
	std::stringstream output;
	output<<"Det"<<det;
	return output.str();
}