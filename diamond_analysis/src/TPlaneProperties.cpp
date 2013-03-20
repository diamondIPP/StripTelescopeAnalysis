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

std::vector<UInt_t> TPlaneProperties::getSiliconPlaneVector() {
	std::vector<UInt_t> vecPlanes;
	for (UInt_t i=0;i<getNSiliconPlanes();i++)
		vecPlanes.push_back(i);
	return vecPlanes;
}

std::string TPlaneProperties::getStringForDetector(int det){
	switch(det){
	case 0: return "D0X";
	case 1: return "D0Y";
	case 2: return "D1X";
	case 3: return "D1Y";
	case 4: return "D2X";
	case 5: return "D2Y";
	case 6: return "D3X";
	case 7: return "D3Y";
	case 8: return "Dia";
	default: return "Invalid";
	}
	return "Invalid";
}
