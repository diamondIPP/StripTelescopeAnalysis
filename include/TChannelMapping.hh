/*
 * TChannelMapping.hh
 *
 *  Created on: Feb 17, 2012
 *      Author: bachmair
 */

#ifndef TCHANNELMAPPING_HH_
#define TCHANNELMAPPING_HH_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include "TPlaneProperties.hh"
#include "TObject.h"

class TChannelMapping:public TObject {
public:
	typedef std::map<UInt_t,UInt_t> channelContainer;
	TChannelMapping(UInt_t nChannels=TPlaneProperties::getNChannelsDiamond());
	TChannelMapping(std::vector<UInt_t> channelMapping);
	TChannelMapping(std::vector<int> channelMapping);
	TChannelMapping(const TChannelMapping &rhs);
	virtual ~TChannelMapping();
	void changeMapping(UInt_t vaChNo,UInt_t detChNo);
	void changeMapping(std::vector<UInt_t> channelMapping);
	UInt_t getDetChannelNo(UInt_t VaChNo);
	UInt_t getVAChannelNo(UInt_t detChNo);
	void PrintMapping();
private:
	void showMappingPair( const  std::pair<UInt_t,UInt_t>  p );
	void showMapping( const  channelContainer::iterator  p );
	channelContainer mapVaToDet;
	std::vector<UInt_t> vecVAChNo;
	std::vector<UInt_t> vecdetChNo;
  ClassDef(TChannelMapping,1);


};

#endif /* TCHANNELMAPPING_HH_ */
