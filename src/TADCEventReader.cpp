/*
 * TADCEvent.cpp
 *
 *  Created on: 01.08.2011
 *      Author: Felix Bachmair
 */

#include "TADCEventReader.hh"
using namespace std;

TADCEventReader::TADCEventReader(string FileName) {
	verbosity=0;
	current_event = 0;
	store_threshold = 2;
	tree =NULL;
	file=NULL;
	/*run_number=RunNumber;
	event_number=EventNumber;*/
	SetTree(FileName);//tree);
	initialiseTree();
	if(!this->isOK()){
		cout<<"TADCEventReader::TADCEventReader is not correctly initialized.. EXIT PROGRAM"<<endl;
		exit(-1);
	}
	if (verbosity) cout<<"tree Entries():"<<tree->GetEntries()<<endl;
	this->GetEvent(0);
	this->fileName=FileName;
}

TADCEventReader::~TADCEventReader() {
	cout<< "deleting instance of TADCEventReader"<<endl;
	delete tree;
	delete file;
}

bool TADCEventReader::SetTree(string fileName){//TTree *tree){
	if(tree!=NULL) tree->Delete();
	if(file!=NULL) file->Delete();
	tree=NULL;
	file=NULL;
	cout<<"load File: \""<<fileName<<"\""<<endl;
	file = new TFile(fileName.c_str());
	file->GetObject("tree",tree);
	if(tree==NULL)
	{
		cout<<"tree does not exist, looking for other TTree objects."<<endl;
		tree=(TTree*)getTreeName();
	}

	if(tree!=NULL){
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
	cout<<"TADCEventReader::isOK \""<<tree<<"\""<<endl;
	return (tree!=NULL);
}
void TADCEventReader::SetBranchAddresses(){
	//Event Header Branches
	tree->SetBranchAddress("RunNumber",&run_number);
	tree->SetBranchAddress("EventNumber",&event_number);
	tree->SetBranchAddress("StoreThreshold",&store_threshold);
//	tree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);
	tree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);

	//Telescope Data Branches
	tree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
	tree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
	tree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
	tree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
	tree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
	tree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
	tree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
	tree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
	tree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
	tree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
	tree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
	tree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
	tree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
	tree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
	tree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
	tree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
	tree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
	tree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
	//tree->SetBranchAddress("Det_ADC",&Det_ADC[0][0]);
	tree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
	tree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
	tree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
	tree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
	tree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
	tree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
	tree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
	tree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
	tree->SetBranchAddress("Dia_ADC",&Dia_ADC);
	tree->SetBranchAddress("DiaADC",&Dia_ADC);
	tree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
	tree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
	tree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
	tree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
	tree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
	tree->SetBranchAddress("D2Y_PedMean",&Det_PedMean[5]);
	tree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
	tree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
	tree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
	tree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
	tree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
	tree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
	tree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
	tree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
	tree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
	tree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
	tree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
	tree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
	tree->SetBranchAddress("PedestalMean",&pedestalMean);
	tree->SetBranchAddress("PedestalSigma",&pedestalSigma);
}

void TADCEventReader::initialiseTree(){
	current_event = 0;
	tree->GetEvent(current_event);
	cout<< "Loaded first event in PedTree: "<<event_number<<endl;
	cout<< "RunNumber is: "<<run_number<<endl;
	cout<< "StoreThreshold is: "<<store_threshold<<endl;
}

bool TADCEventReader::GetNextEvent(){
	if(tree==NULL) return true;
	if(current_event+1<tree->GetEntries()){
		current_event++;
		tree->GetEvent(current_event);
		if(verbosity>=2)
			cout<<"Got next Event, new event number: "<<current_event<<endl;
		return true;
	}
	return false;
}
bool TADCEventReader::GetEvent(UInt_t EventNumber){
	if(tree==NULL) return false;
	if(EventNumber<tree->GetEntries()){
			current_event=EventNumber;
			tree->GetEvent(current_event);
			if(verbosity>=2)
				cout<<"Got Event: "<<current_event<<endl;
			return true;
		}
	return false;
}

Long64_t TADCEventReader::GetEntries(){
	if (verbosity>=2) {
		cout<<"TADCEventReader::GetEntries:"<<tree<<flush;
		cout<<" "<<tree->GetEntries();
	}
	if(tree!=NULL)
		return tree->GetEntries();
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

UChar_t TADCEventReader::getDet_ADC(UInt_t det , UInt_t ch) const
{
    return Det_ADC[det][ch];
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

UShort_t TADCEventReader::getDia_ADC(UInt_t ch) const
{
    return Dia_ADC[ch];
}

UInt_t TADCEventReader::getEvent_number() const
{
    return event_number;
}

TTree * TADCEventReader::getPedTree() const
{
    return tree;
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

int TADCEventReader::hasTree(){
	if(file==NULL)return -1;
	  TIter nextkey(file->GetListOfKeys());
	  TKey *key;
	  int hasATree = 0;
	  while ((key = (TKey*)nextkey()))
	    {
	      TObject *obj = key->ReadObj();
	      if ((obj->IsA()->InheritsFrom("TTree"))){
	    	  hasATree++;}
	    }
	  return hasATree;
}

std::string TADCEventReader::getStringForPlane(int i)
{
	switch(i){
	case 0: return "D0X";
	case 1: return "D0Y";
	case 2: return "D1X";
	case 3: return "D1Y";
	case 4: return "D2X";
	case 5: return "D2Y";
	case 6: return "D3X";
	case 7: return "D3Y";
	case 8: return "Dia";
	}
}

TFile *TADCEventReader::getFile() const
{
	return this->file;
}

Float_t TADCEventReader::getPedestalMean(UInt_t det, UInt_t ch)
{
	return this->pedestalMean[det][ch];
}

Float_t TADCEventReader::getPedestalSigma(UInt_t det, UInt_t ch)
{
	return this->pedestalSigma[det][ch];
}

TObject* TADCEventReader::getTreeName(){
	cout<<"TADCEventReader::getTreeName:"<<endl;
	if(file==NULL) exit(-1);
	TObject *obj=NULL;
	int hastree=hasTree();
		cout<< "File has "<<hastree<<" trees"<<endl;
//	if(hastree!=1)
//		return obj;
	TIter nextkey(file->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)nextkey()))
	{
		obj = key->ReadObj();
		if ((obj->IsA()->InheritsFrom("TTree"))){
			return obj;
		}
		else
			obj=NULL;
	}
	return NULL;
}

TTree *TADCEventReader::getTree() const
{
    return tree;
}

std::string TADCEventReader::getFilePath(){
	return this->fileName;
}



