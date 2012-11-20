/*
 * TChannelMapping.cpp
 *
 *  Created on: Feb 17, 2012
 *      Author: bachmair
 */

#include "../include/TChannelMapping.hh"

ClassImp(TChannelMapping);
using namespace std;
TChannelMapping::TChannelMapping(UInt_t nChannels)
{
	cout<<"ChannelMapping:Konstruktor nChannels "<<nChannels<<flush;
	for(UInt_t ch =0;ch<nChannels;ch++)
		changeMapping(ch,ch);
	cout<<"PRINT:"<<flush;
//	this->PrintMapping();
	cout<<"DONE"<<endl;
}

TChannelMapping::TChannelMapping(std::vector<UInt_t> channelMapping)
{
	cout<<"TChannelMapping:: Constructor vect UInt..."<<flush;
	for(UInt_t ch=0;ch<channelMapping.size();ch++){
	}
//		this->mapVaToDet[ch]=channelMapping.at(ch);
//	this->PrintMapping();
	cout<<"DONE"<<endl;
}
TChannelMapping::TChannelMapping(std::vector<int> channelMapping)
{
	cout<<"TChannelMapping:: Constructor vect Int..."<<channelMapping.size()<<flush;
	for(UInt_t ch=0;ch<channelMapping.size();ch++)
		changeMapping(ch,channelMapping.at(ch));
//	this->PrintMapping();
	cout<<"DONE"<<endl;
}

TChannelMapping::TChannelMapping(const TChannelMapping & rhs)
{
	cout<<"TChannelMapping:: COPY Constructor..."<<rhs.mapVaToDet.size()<<flush;
	std::map<UInt_t,UInt_t>::const_iterator it;
	for ( it=rhs.mapVaToDet.begin(); it != rhs.mapVaToDet.end(); it++){
		changeMapping((*it).first,(*it).second);
	}
	this->PrintMapping();
	cout<<"DONE"<<endl;

}

TChannelMapping::~TChannelMapping() {
	// TODO Auto-generated destructor stub
}

void TChannelMapping::changeMapping(UInt_t vaChNo, UInt_t detChNo)
{
	if(mapVaToDet.find(vaChNo)!=mapVaToDet.end())
		mapVaToDet.at(vaChNo)=detChNo;
	else
		this->mapVaToDet.insert(std::pair<UInt_t,UInt_t>(vaChNo,detChNo));//[ch]=ch;
//	mapVaToDet[vaChNo]=detChNo;
}


void TChannelMapping::changeMapping(vector<UInt_t> channelMapping){
	for(UInt_t ch=0;ch<channelMapping.size();ch++);
//		this->mapVaToDet[ch]=channelMapping.at(ch);
//	this->PrintMapping();
}


UInt_t TChannelMapping::getDetChannelNo(UInt_t vaChNo)
{
	channelContainer::iterator it;
	for ( it=mapVaToDet.begin() ; it != mapVaToDet.end(); it++ ){
//		cout<<"\t"<<setw(3)<<(*it).first<<"-->"<<setw(3)<<(*it).second<<"\n";
//		showMapping(it);
		if((*it).first==vaChNo)break;
	}
//	cout<<"###"<<flush;
	if(it==mapVaToDet.end()){
		cerr<<vaChNo<<" does not exist:"<<endl;
		return vaChNo;

	}
	return (*it).second;
//
//	return mapVaToDet[vaChNo];
	return 999;
}



UInt_t TChannelMapping::getVAChannelNo(UInt_t detChNo)
{
	std::map<UInt_t,UInt_t>::iterator it;
	for ( it=mapVaToDet.begin() ; it != mapVaToDet.end(); it++ )
		if((*it).second==detChNo)break;
	if(it==mapVaToDet.end()){
		cerr<<detChNo<<" does not exist:"<<endl;
		return -1;
	}
	return (*it).first;
	return 999;
}




void TChannelMapping::PrintMapping()
{
//	std::for_each( mapVaToDet.begin(), mapVaToDet.end(), showMappingPair);
	cout<< "TChannelMapping::PrintMapping: Channel Map: of "<<mapVaToDet.size()<<" Channels"<<endl;
	channelContainer::iterator it;
	for ( it=mapVaToDet.begin() ; it != mapVaToDet.end(); it++ ){
		UInt_t ch = (*it).first;
		showMapping(it);
//		cout<<"\t"<<setw(3)<<(*it).first<<"-->"<<
//		getDetChannelNo(ch)<<
////		getDetChannelNo( (it*).first )<<
//		setw(3)<<(*it).second<<"\n";
	}
	cout<<endl;
}

void TChannelMapping::showMappingPair( const std::pair<UInt_t,UInt_t> p )
{
	try{
		std::cout	  << p.first << "  -> " <<getDetChannelNo(p.first) << std::endl;
	}
	catch( std::exception &Fehler )
	{
//		std::cerr << "Fehler: <" << Fehler.what() << "> in showMapping(pair)\n"<<flush;
	}
	catch( ... )
	{
		std::cerr << "unbekannter Fehler showMapping(pair)\n"<<flush;
}
}

void TChannelMapping::showMapping(const channelContainer::iterator p)
{
	try{
	  std::cout	  << (*p).first << "  -> " << getDetChannelNo((*p).first)  << std::endl;
	}
	catch( std::exception &Fehler )
	{
		std::cerr << "Fehler: <" << Fehler.what() << "> in showMapping(iterator)\n"<<flush;
	}
	catch( ... )
	{
		std::cerr << "unbekannter Fehler showMapping(iterator)\n"<<flush;
	}

}




