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
	verbosity = 0;
}

TDiamondPattern::~TDiamondPattern() {
	// TODO Auto-generated destructor stub
}

void TDiamondPattern::loadStandardPitchWidthSettings() {
	loadPitchWidthSettings(TPlaneProperties::getStripDistance());
	bLoadedStandardPitchWidthSettings=true;
}


void TDiamondPattern::loadPitchWidthSettings(Float_t pitchWidth) {
	resetPattern();
	addPattern(pitchWidth,0,0,channelToMetricConversion.size()-1);
	standardPW = pitchWidth;
}

void TDiamondPattern::resetPattern() {
	initialiseVector();
	beginOfInterval.clear();
	endOfInterval.clear();
	firstChannelOfInterval.clear();
	nChannelsOfInterval.clear();
	pitchWidth.clear();
	bLoadedStandardPitchWidthSettings=false;
}

bool TDiamondPattern::addPattern(Float_t pitchWidth, Float_t startPosition, UInt_t firstChannel, UInt_t lastChannel) {
	if(firstChannel<0||lastChannel>=channelToMetricConversion.size()||(int)lastChannel-(int)firstChannel<0){
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
	Float_t endPos = channelToMetricConversion[lastChannel];
	endOfInterval.push_back(endPos);
	firstChannelOfInterval.push_back(firstChannel);
	UInt_t deltaChannels = lastChannel-firstChannel;

	nChannelsOfInterval.push_back(deltaChannels);
	this->pitchWidth.push_back(pitchWidth);
	cout<<"added: "<<getInterval(beginOfInterval.size()-1).first<<"-"<<getInterval(beginOfInterval.size()-1).second<<endl;
	return retVal;
}

bool TDiamondPattern::isValidChannelPosition(Float_t channel) {
	for ( UInt_t i = 0; i <nChannelsOfInterval.size(); i++){
		if (firstChannelOfInterval[i]-.5 < channel && channel < firstChannelOfInterval[i]+nChannelsOfInterval[i] +.5)
			return true;
	}
	return false;
}

bool TDiamondPattern::isValidCluster(TCluster* cluster) {
	for (UInt_t cl = 0; cl < cluster->size(); cl ++)
		if (cluster->isHit(cl)||cluster->isSeed(cl))
			if ( !isValidChannelPosition( cluster->getChannel(cl) ) )
				return false;
	return true;
}

void TDiamondPattern::initialiseVector() {
	for (UInt_t i=0;i<channelToMetricConversion.size();i++)
		channelToMetricConversion[i] = N_INVALID;
}


void TDiamondPattern::Print() {
	cout<<"Diamond Detector Pattern \n"<<endl;
	//	cout<<"Valid: "<< this->hasInvalidIntervals()<<endl;
	cout<<"  ch | Pos [um]"<<endl;
	cout<<"-----+---------"<<endl;
	for(UInt_t i = 0; i<channelToMetricConversion.size();i++){
		cout<<" "<<setw(3)<<i<<" | "<< setw(6)<<convertChannelToMetric(i)<<endl;
	}
	cout<<endl;

	cout<<"nChannelsOfInterval:   "<<nChannelsOfInterval.size()<<endl;
	cout<<"beginOfInterval:       "<<beginOfInterval.size()<<endl;
	cout<<"endOfInterval:         "<<endOfInterval.size()<<endl;
	cout<<"firstChannelOfInterval "<<firstChannelOfInterval.size()<<endl;
	cout<<"pitchWidth:            "<<pitchWidth.size()<<endl;
	//	char t;
	//	cin>>t;
}


UInt_t TDiamondPattern::getNIntervals(){
	if (hasInvalidIntervals())
		return 0;
	return nChannelsOfInterval.size();
}

bool TDiamondPattern::hasInvalidIntervals(){
	bool retVal;
	retVal = nChannelsOfInterval.size() == beginOfInterval.size();
	retVal = retVal && (beginOfInterval.size() == endOfInterval.size());
	retVal = retVal && (endOfInterval.size() == firstChannelOfInterval.size());
	if(!retVal){
		cerr << " ERROR!! *********************************"<<endl;
		cout <<"hasInvalidIntervals! "<< nChannelsOfInterval.size() << " " << beginOfInterval.size() << " " << endOfInterval.size() << " " << firstChannelOfInterval.size() <<endl;
		//		this->Print();
		exit(-1);
	}

	return !retVal;
}
void TDiamondPattern::showPatterns(){
	cout<< "There are "<<getNIntervals()<<"/"<<getNPatterns()<<"  Patterns:"<<endl;
	for (UInt_t i=0;i<getNIntervals();i++){
		std::pair<Int_t, Int_t> interval = getInterval(i);
		cout<<"Pattern "<<i<<":\t"<<interval.first<<"-"<<interval.second<<endl;
	}
}

Float_t TDiamondPattern::convertChannelToMetric(Float_t channel) {
	if(channel>=channelToMetricConversion.size()||channel<0){
		//		cout<<" TDiamondPattern::convertChannelToMetric "<<channel<<" INVALID"<<endl;
		return N_INVALID;
	}
	int leftChannel = (UInt_t)channel;
	int rightChannel = leftChannel+1;
	Float_t leftPosition =  getChannelToMetric(leftChannel);
	Float_t rightPosition = getChannelToMetric(rightChannel);
	Float_t position = leftPosition + (rightPosition-leftPosition)*(channel-leftChannel);
	if(leftChannel==channel)
		return leftPosition;
	if (leftPosition == N_INVALID || rightPosition == N_INVALID){
		if(verbosity) cerr<<"One of the channels is invalid: "<<leftChannel<<"-"<<rightChannel<<"\t"<<leftPosition<<"-"<<rightPosition<<endl;
		return N_INVALID;
	}
	return position;
}


Int_t TDiamondPattern::getPatternOfHit(Float_t metric) {
	if(hasInvalidIntervals()){

		if(verbosity)  cout<<" TDiamondPattern::convertMetricToChannel "<<metric<<" INVALID INTERVALS"<<endl;
		return N_INVALID;
	}
	for(UInt_t i = 0; i < getNIntervals();i++){
		if(i<beginOfInterval.size() && i < endOfInterval.size()){
			if( beginOfInterval[i] <= metric && metric <= endOfInterval[i] )
				return i;
		}
	}
	return N_INVALID;
}

Float_t TDiamondPattern::convertMetricToChannel(Float_t metric) {
	if(verbosity>5)
		cout<<"[TDiamondPattern::convertMetricToChannel]"<<metric<<endl;
	if(hasInvalidIntervals()){
		if(verbosity)
			cout<<" TDiamondPattern::convertMetricToChannel "<<metric<<" INVALID INTERVALS"<<endl;
		return N_INVALID;
	}
	for(UInt_t i = 0; i < getNIntervals();i++){
		if(i<beginOfInterval.size() && i < endOfInterval.size() && i < pitchWidth.size()){
			Float_t pw = pitchWidth[i];
			Float_t begin = beginOfInterval[i];
			Float_t end = endOfInterval[i];
			begin -= pw/2;
			end += pw/2;

//			if(metric >= 3660 && metric <= 3715)
//				cout<<metric <<"look at" <<begin<<"-"<<end<<" "<<pw<<endl;
			if( begin <= metric && metric <= end ){
				Float_t retval =convertMetricToChannel(metric,i);
//				if(metric >= 3660 && metric <= 3715)
//					cout<<"convert"<<i<<":\t"<<retval<<endl;
				return retval;
			}
		}
		//		else return N_INVALID;
	}
	return N_INVALID;
}


Float_t TDiamondPattern::convertMetricToRelativeMetric(Float_t metric,UInt_t interval){
    if (interval>=getNIntervals())
        return N_INVALID;
    Float_t leftPos,rightPos;
    Int_t firstChannel = firstChannelOfInterval[interval]-1;
    Int_t lastChannel = firstChannelOfInterval[interval] + nChannelsOfInterval[interval];
//  if(metric >= 3660 && metric <= 3715) cout<<"[TDiamondPattern::convertMetricToChannel]"<<metric<<" in " << interval <<": "<<firstChannel<<"-"<<lastChannel<<endl;
    for(Int_t ch = firstChannel; ch <= lastChannel;ch++){
//      if(metric >= 3660 && metric <= 3715)
//          cout<<"check channel "<<ch<<" "<<flush;
        if(ch == firstChannel)
            leftPos = getChannelToMetric(ch+1) - pitchWidth[interval];
        else
            leftPos = getChannelToMetric(ch);
        if (ch == lastChannel)
            rightPos = getChannelToMetric(ch) + pitchWidth[interval];
        else
            rightPos = getChannelToMetric(ch+1);
//      if(metric >= 3660 && metric <= 3715)
//          cout<<" "<<leftPos<<"-"<<rightPos<<endl;
        if ( leftPos <= metric && metric <= rightPos ){
            if(metric-leftPos < rightPos-metric)
                return metric - leftPos;
            else
                return metric - rightPos;
        }
    }
    return N_INVALID;

}
Float_t TDiamondPattern::convertMetricToRelativeMetric(Float_t metric) {
    if(verbosity>5)
        cout<<"[TDiamondPattern::convertMetricToRelativeMetric]"<<metric<<endl;
    if(hasInvalidIntervals()){
        if(verbosity)
            cout<<" TDiamondPattern::convertMetricToRelativeMetric "<<metric<<" INVALID INTERVALS"<<endl;
        return N_INVALID;
    }
    for(UInt_t i = 0; i < getNIntervals();i++){
        if(i<beginOfInterval.size() && i < endOfInterval.size() && i < pitchWidth.size()){
            Float_t pw = pitchWidth[i];
            Float_t begin = beginOfInterval[i];
            Float_t end = endOfInterval[i];
            begin -= pw/2;
            end += pw/2;

            if( begin <= metric && metric <= end ){
                Float_t retval = convertMetricToRelativeMetric(metric,i);
                return retval;
            }
        }
        //      else return N_INVALID;
    }
    return N_INVALID;
}

Float_t TDiamondPattern::convertMetricToChannel(Float_t metric,UInt_t interval) {
	if (interval>=getNIntervals())
		return N_INVALID;
	Float_t leftPos,rightPos;
	Int_t firstChannel = firstChannelOfInterval[interval]-1;
	Int_t lastChannel = firstChannelOfInterval[interval] + nChannelsOfInterval[interval];
//	if(metric >= 3660 && metric <= 3715) cout<<"[TDiamondPattern::convertMetricToChannel]"<<metric<<" in " << interval <<": "<<firstChannel<<"-"<<lastChannel<<endl;
	for(Int_t ch = firstChannel; ch <= lastChannel;ch++){
//		if(metric >= 3660 && metric <= 3715)
//			cout<<"check channel "<<ch<<" "<<flush;
		if(ch == firstChannel)
			leftPos = getChannelToMetric(ch+1) - pitchWidth[interval];
		else
			leftPos = getChannelToMetric(ch);
		if (ch == lastChannel)
			rightPos = getChannelToMetric(ch) + pitchWidth[interval];
		else
			rightPos = getChannelToMetric(ch+1);
//		if(metric >= 3660 && metric <= 3715)
//			cout<<" "<<leftPos<<"-"<<rightPos<<endl;
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


std::pair<Int_t,Int_t> TDiamondPattern::getInterval(UInt_t pattern){
    Int_t first = -1;
    Int_t last = -1;
    if(pattern >= 0 && pattern<getNIntervals()){
        first= firstChannelOfInterval.at(pattern);
        last = nChannelsOfInterval.at(pattern) + first;
    }
    return make_pair(first,last);
}
std::pair<Int_t, Int_t> TDiamondPattern::getTotalInterval(){

    Int_t first = TPlaneProperties::getNChannelsDiamond();
    Int_t last = 0;
    for (UInt_t area = 0; area < getNIntervals();area ++){
        pair<Int_t,Int_t> pattern  = getInterval(area);
        if(pattern.first < first) first = pattern.first;
        if(pattern.second > last) last = pattern.second;
    }

    return make_pair(first,last);
}


std::pair<int,int> TDiamondPattern::getPatternChannels(UInt_t pattern){
	int firstChannel=-1;
	int lastChannel=-1;
	pattern--;
	if(pattern<firstChannelOfInterval.size()&&pattern<nChannelsOfInterval.size()){
		firstChannel = firstChannelOfInterval[pattern];
		lastChannel = nChannelsOfInterval[pattern];
		lastChannel+=firstChannel;
	}
	return make_pair(firstChannel,lastChannel);
}

Int_t TDiamondPattern::getClusterPattern(TCluster* cluster) {
    Int_t first = cluster->getFirstHitChannel();
    Int_t last = cluster->getLastHitChannel();
    if(verbosity>7)cout<<"first: "<<first<<" - last: "<<last<<endl;
    for(UInt_t i = 0; i< getNPatterns();i++){
        pair<int,int> channels = getPatternChannels(i+1);
        if(verbosity>7)cout<<channels.first<<"<="<<first<<" "<<last<<"<="<<channels.second<<endl;
        if(channels.first<=first&&last<=channels.second)
            return i+1;
    }
    return -1;
}

Int_t TDiamondPattern::getPatternOfChannel(Int_t ch) {
    for(UInt_t i = 0; i< getNPatterns();i++){
        pair<int,int> channels = getPatternChannels(i+1);
        if(verbosity>7)cout<<channels.first<<"<="<<ch<<"<="<<channels.second<<endl;
        if(channels.first<=ch&&ch<=channels.second)
            return i+1;
    }
    return -1;
}
