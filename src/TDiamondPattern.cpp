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
}

TDiamondPattern::~TDiamondPattern() {
	// TODO Auto-generated destructor stub
}

void TDiamondPattern::loadStandardPitchWidthSettings() {
	loadPitchWidthSettings(TPlaneProperties::getStripDistance());
}

void TDiamondPattern::loadPitchWidthSettings(Float_t pitchWidth) {
	resetPattern();
	addPattern(pitchWidth,0,0,channelToMetricConversion.size()-1);
}

void TDiamondPattern::resetPattern() {
	initialiseVector();
}

bool TDiamondPattern::addPattern(Float_t pitchWidth, Float_t startPosition, UInt_t firstChannel, UInt_t lastChannel) {
	if(firstChannel<0||lastChannel>=channelToMetricConversion.size()){
		cout<< "Want to create pattern with invalid channel No: "<<firstChannel<<"-"<<lastChannel<<endl;
		return false;
	}
	cout<< "Adding new Pattern with a pitchWidth of "<<pitchWidth<<" um @ "<<startPosition << " um, Channels: "<<firstChannel <<" - "<<lastChannel<<endl;
	bool retVal = true;
	for(UInt_t i = firstChannel;i<=lastChannel;i++){
		if(channelToMetricConversion[i]!=0) retVal = false;
		channelToMetricConversion[i] = startPosition + pitchWidth*(i-lastChannel);
	}
	beginOfInterval.push_back(startPosition);
	endOfInterval.push_back(channelToMetricConversion[lastChannel]);
	firstChannelOfInterval.push_back(firstChannel);
	return retVal;
}


void TDiamondPattern::initialiseVector() {
	for (UInt_t i=0;i<channelToMetricConversion.size();i++)
			channelToMetricConversion[i] = 0;
}


void TDiamondPattern::Print() {
	cout<<"Diamond Detector Pattern \n"<<endl;
	cout<<"  ch | Pos [um]"<<endl;
	cout<<"-----+---------"<<endl;
	for(UInt_t i = 0; i<channelToMetricConversion.size();i++){
		cout<<" "<<setw(3)<<i<<" | "<< setw(6)<<channelToMetricConversion[i]<<endl;
	}
	cout<<endl;
	char t;
	cin>>t;
}

Float_t TDiamondPattern::convertChannelToMetric(Float_t channel) {
	return -1;
}

Float_t TDiamondPattern::convertMetricToChannel(Float_t metric) {
	return -1;
}
