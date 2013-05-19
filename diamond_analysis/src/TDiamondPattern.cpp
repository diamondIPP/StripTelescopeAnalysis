/**
 * @file TDiamondPattern.cpp
 *
 * @date Jan 9, 2013
 * @author bachmair
 * @description
 */

#include "../include/TDiamondPattern.hh"
using namespace std;
ClassImp(TDiamondPattern);
TDiamondPattern::TDiamondPattern() {
	// TODO Auto-generated constructor stub
	channelToMetricConversion.resize(TPlaneProperties::getNChannelsDiamond());
	bLoadedStandardPitchWidthSettings=false;
	standardPW = TPlaneProperties::getStripDistance();
}

TDiamondPattern::~TDiamondPattern() {
	// TODO Auto-generated destructor stub
}

void TDiamondPattern::loadStandardPitchWidthSettings() {
	loadPitchWidthSettings(TPlaneProperties::getStripDistance());
	bLoadedStandardPitchWidthSettings=true;
}

std::pair<Int_t,Int_t> TDiamondPattern::getInterval(UInt_t pattern){
	Int_t first = -1;
	Int_t last = -1;
	if(pattern >= 0 && pattern<getNIntervals()){
		first= firstChannelOfInterval.at(pattern);
		last = nChannelsOfInterval.at(pattern) + first -1;
	}
	return make_pair(first,last);
}

void TDiamondPattern::loadPitchWidthSettings(Float_t pitchWidth) {
	resetPattern();
	addPattern(pitchWidth,0,0,channelToMetricConversion.size()-1);
	standardPW =pitchWidth;
}

void TDiamondPattern::resetPattern() {
	initialiseVector();
	beginOfInterval.clear();
	endOfInterval.clear();
	firstChannelOfInterval.clear();
	nChannelsOfInterval.clear();
	bLoadedStandardPitchWidthSettings=false;
}

bool TDiamondPattern::addPattern(Float_t pitchWidth, Float_t startPosition, UInt_t firstChannel, UInt_t lastChannel) {
	if(firstChannel<0||lastChannel>=channelToMetricConversion.size()||lastChannel-firstChannel<0){
		cout<< "Want to create pattern with invalid channel No: "<<firstChannel<<"-"<<lastChannel<<endl;
		return false;
	}
	bool retVal = true;
	Float_t endPosition = startPosition+pitchWidth*(lastChannel-firstChannel);
	//check If interval overlaps with existing one
	UInt_t i;
	for( i=0;i<nChannelsOfInterval.size() && retVal == true;i++){
		if(beginOfInterval[i]<startPosition && startPosition<endOfInterval[i])
			retVal = false;
		if(beginOfInterval[i]<endPosition && endPosition<endOfInterval[i])
			retVal = false;
		if(firstChannelOfInterval[i] <= firstChannel && firstChannel <= firstChannelOfInterval[i] + nChannelsOfInterval[i])
			retVal = false;
		if(firstChannelOfInterval[i] <= lastChannel && lastChannel <= firstChannelOfInterval[i] + nChannelsOfInterval[i])
			retVal = false;
	}
	if (retVal == false){
		cout<<"Couldn't create Pattern, pattern overlaps with another one."<<i<<endl;
//		this->Print();
		return retVal;
	}
	cout<< "Adding new Pattern with a pitchWidth of "<<pitchWidth<<" um @ "<<startPosition << " um, Channels: "<<firstChannel <<" - "<<lastChannel<<endl;
	Float_t pos = startPosition;
	for(UInt_t i = firstChannel;i<=lastChannel;i++){
		if(channelToMetricConversion[i] != N_INVALID) retVal = false;
		//		cout<<"add "<<i<<": "<<pos<<endl;
		channelToMetricConversion[i] = pos;
		pos +=pitchWidth;
	}
	beginOfInterval.push_back(startPosition);
	endOfInterval.push_back(channelToMetricConversion[lastChannel]);
	firstChannelOfInterval.push_back(firstChannel);
	nChannelsOfInterval.push_back(lastChannel-firstChannel);
	this->pitchWidth.push_back(pitchWidth);
	return retVal;
}


void TDiamondPattern::initialiseVector() {
	for (UInt_t i=0;i<channelToMetricConversion.size();i++)
		channelToMetricConversion[i] = N_INVALID;
}


void TDiamondPattern::Print() {
	cout<<"Diamond Detector Pattern \n"<<endl;
	cout<<"  ch | Pos [um]"<<endl;
	cout<<"-----+---------"<<endl;
	for(UInt_t i = 0; i<channelToMetricConversion.size();i++){
		cout<<" "<<setw(3)<<i<<" | "<< setw(6)<<convertChannelToMetric(i)<<endl;
	}
	cout<<endl;
	//	char t;
	//	cin>>t;
}

void TDiamondPattern::showPatterns(){
	cout<< "There are "<<getNIntervals()<<"/"<<getNPatterns()<<"  Patterns:"<<endl;
	for (UInt_t i=0;i<getNIntervals();i++){
		std::pair<Int_t, Int_t> interval = getInterval(i);
		cout<<"Pattern "<<i<<":\t"<<interval.first<<"-"<<interval.second<<endl;
	}
}

Float_t TDiamondPattern::convertChannelToMetric(Float_t channel) {
	if(channel>=channelToMetricConversion.size()||channel<0)
		return N_INVALID;
	int leftChannel = (UInt_t)channel;
	int rightChannel = leftChannel+1;
	Float_t leftPosition =  getChannelToMetric(leftChannel);
	Float_t rightPosition = getChannelToMetric(rightChannel);
	Float_t position = leftPosition + (rightPosition-leftPosition)*(channel-leftChannel);
	if(leftChannel==channel)
		return leftPosition;
	if (leftPosition == N_INVALID || rightPosition == N_INVALID){
		//cerr<<"One of the channels is invalid: "<<leftChannel<<"-"<<rightChannel<<"\t"<<leftPosition<<"-"<<rightPosition<<endl;
		return N_INVALID;
	}
	return position;
}

Float_t TDiamondPattern::convertMetricToChannel(Float_t metric) {
	for(UInt_t i = 0; i < getNIntervals();i++){
		if( beginOfInterval[i] <= metric && metric <= endOfInterval[i] )
			return convertMetricToChannel(metric,i);
	}
	return N_INVALID;
}

Float_t TDiamondPattern::convertMetricToChannel(Float_t metric,UInt_t interval) {
	if (interval>=getNIntervals())
		return N_INVALID;
	for(UInt_t ch=firstChannelOfInterval[interval];ch < firstChannelOfInterval[interval] + nChannelsOfInterval[interval];ch++){
		Float_t leftPos = getChannelToMetric(ch);
		Float_t rightPos = getChannelToMetric(ch+1);
		if ( leftPos <= metric && metric <= rightPos )
			return (metric-leftPos)/(rightPos-leftPos)+ch;
	}
	return N_INVALID;
}

Float_t TDiamondPattern::getChannelToMetric(UInt_t ch){
	if(ch>=0&&ch<channelToMetricConversion.size())
		return channelToMetricConversion[ch];
	return N_INVALID;
}

Float_t TDiamondPattern::getPitchWidth(UInt_t area){
	if(area<beginOfInterval.size())
		return this->pitchWidth.at(area);
	else
		return this->standardPW;
}
