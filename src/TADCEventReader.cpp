/*
 * TADCEvent.cpp
 *
 *  Created on: 01.08.2011
 *      Author: Felix Bachmair
 */

#include "TADCEventReader.hh"
using namespace std;

TADCEventReader::TADCEventReader(string fileName) {
	verbosity=0;
	current_event = 0;
	store_threshold = 2;
	PedTree =NULL;
	PedFile=NULL;
	/*run_number=RunNumber;
	event_number=EventNumber;*/
	SetTree(fileName);//tree);
	initialiseTree();
	if(!this->isOK()){
		cout<<"TADCEventReader::TADCEventReader is not correctly initialized.. EXIT PROGRAM"<<endl;
		exit(-1);
	}
	if (verbosity) cout<<"PedTree Entries():"<<PedTree->GetEntries()<<endl;
	this->GetEvent(0);
}

TADCEventReader::~TADCEventReader() {
	cout<< "deleting instance of TADCEventReader"<<endl;
	delete PedTree;
	delete PedFile;
}

bool TADCEventReader::SetTree(string fileName){//TTree *tree){
	if(PedTree!=NULL) PedTree->Delete();
	if(PedFile!=NULL) PedFile->Delete();
	PedTree=NULL;
	PedFile=NULL;
	PedFile = new TFile(fileName.c_str());
	PedFile->GetObject("PedTree",PedTree);

	if(PedTree!=NULL){
		SetBranchAddresses();
		return true;
	}
	else {
		if(verbosity)
			cerr<<"The given Tree is not correct, please check InputTree"<<endl;
		return false;
	}
}
bool TADCEventReader::isOK(){
	cout<<"TADCEventReader::isOK \""<<PedTree<<"\""<<endl;
	return (PedTree!=NULL);
}
void TADCEventReader::SetBranchAddresses(){
	//Event Header Branches
	PedTree->SetBranchAddress("RunNumber",&run_number);
	PedTree->SetBranchAddress("EventNumber",&event_number);
	PedTree->SetBranchAddress("StoreThreshold",&store_threshold);
//	PedTree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);
	PedTree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);

	//Telescope Data Branches
	PedTree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
	PedTree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
	PedTree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
	PedTree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
	PedTree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
	PedTree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
	PedTree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
	PedTree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
	PedTree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
	PedTree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
	PedTree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
	PedTree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
	PedTree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
	PedTree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
	PedTree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
	PedTree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
	PedTree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
	PedTree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
	PedTree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
	PedTree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
	PedTree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
	PedTree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
	PedTree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
	PedTree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
	PedTree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
	PedTree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
	PedTree->SetBranchAddress("Dia_ADC",&Dia_ADC);
	PedTree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
	PedTree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
	PedTree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
	PedTree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
	PedTree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
	PedTree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
	PedTree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
	PedTree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
	PedTree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
	PedTree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
	PedTree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
	PedTree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
	PedTree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
	PedTree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
	PedTree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
	PedTree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
	PedTree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
	PedTree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
}

void TADCEventReader::initialiseTree(){
	current_event = 0;
	PedTree->GetEvent(current_event);
	cout<< "Loaded first event in PedTree: "<<event_number<<endl;
	cout<< "RunNumber is: "<<run_number<<endl;
	cout<< "StoreThreshold is: "<<store_threshold<<endl;
}

bool TADCEventReader::GetNextEvent(){
	if(PedTree==NULL) return true;
	if(current_event+1<PedTree->GetEntries()){
		current_event++;
		PedTree->GetEvent(current_event);
		if(verbosity>=2)
			cout<<"Got next Event, new event number: "<<current_event<<endl;
		return true;
	}
	return false;
}
bool TADCEventReader::GetEvent(UInt_t EventNumber){
	if(PedTree==NULL) return false;
	if(EventNumber<PedTree->GetEntries()){
			current_event=EventNumber;
			PedTree->GetEvent(current_event);
			if(verbosity>=2)
				cout<<"Got Event: "<<current_event<<endl;
			return true;
		}
	return false;
}

Long64_t TADCEventReader::GetEntries(){
	if (verbosity>=2) {
		cout<<"TADCEventReader::GetEntries:"<<PedTree<<flush;
		cout<<" "<<PedTree->GetEntries();
	}
	if(PedTree!=NULL)
		return PedTree->GetEntries();
	else return 0;
}

/*bool TADCEventReader::getCMNEvent_flag() const
{
    return CMNEvent_flag;
}*/

UInt_t TADCEventReader::getCurrent_event() const
{
    return current_event;
}

UChar_t TADCEventReader::getDet_ADC(UInt_t i , UInt_t j) const
{
    return Det_ADC[i][j];
}

UChar_t TADCEventReader::getDet_Channels(UInt_t i , UInt_t j) const
{
	if (i>8||j>255){
		cout<<"TADCEventReader::getDet_Channels not Valid "<<i<<" "<<j<<endl;
		exit (-1);
	}
    return Det_Channels[i][j];
}

UInt_t TADCEventReader::getDet_NChannels(UInt_t i)  const
{
    return Det_NChannels[i];
}

Float_t TADCEventReader::getDet_PedMean(UInt_t i, UInt_t j) const
{
    return Det_PedMean[i][j];
}

Float_t TADCEventReader::getDet_PedWidth(UInt_t i, UInt_t j) const
{
    return Det_PedWidth[i][j];
}

UShort_t TADCEventReader::getDia_ADC(UInt_t i) const
{
    return Dia_ADC[i];
}

UInt_t TADCEventReader::getEvent_number() const
{
    return event_number;
}

TTree * TADCEventReader::getPedTree() const
{
    return PedTree;
}

UInt_t  TADCEventReader::getRun_number() const
{
    return run_number;
}

Float_t  TADCEventReader::getStore_threshold() const
{
    return store_threshold;
}

UInt_t  TADCEventReader::getVerbosity() const
{
    return verbosity;
}

bool  TADCEventReader::getZeroDivisorEvent_flag() const
{
    return ZeroDivisorEvent_flag;
}

