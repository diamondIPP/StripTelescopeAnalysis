/*
 * TADCEvent.cpp
 *
 *  Created on: 01.08.2011
 *      Author: Felix Bachmair
 */

#include "TADCEventReader.hh"
using namespace std;

TADCEventReader::TADCEventReader(string FileName,UInt_t runNumber) {

  cout<<"new TADCEventReader: \n\tpathName:"<<FileName<<"\n\tRunNumber: "<<runNumber<<endl;
	verbosity=0;
  pEvent=NULL;//new TEvent();
	current_event = 0;
	store_threshold = 2;
	tree =NULL;
	file=NULL;
	sys = gSystem;
	if(verbosity>1)cout<<"SYS: "<<sys->pwd()<<endl;
	if(verbosity>1)cout<<"OPEN: "<<FileName<<endl;
	SetTree(FileName);//tree);
	initialiseTree();
	if(!this->isOK()){
		cout<<"TADCEventReader::TADCEventReader is not correctly initialized.. EXIT PROGRAM"<<endl;
		exit(-1);
	}
	if (verbosity) cout<<"tree Entries():"<<tree->GetEntries()<<endl;
	this->LoadEvent(0);
	this->fileName=FileName;

  for(UInt_t det=0;det<9;det++)hEtaIntegral[det]=0;
	LoadEtaDistributions(runNumber);
	pEvent=NULL;//new TEvent();
}

TADCEventReader::~TADCEventReader() {
	cout<< "deleting instance of TADCEventReader"<<flush;
	//delete tree;
	if(file!=0)delete file;
	cout<< "DONE"<<flush;
//	gSystem->Clear();
//	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++)
//	  if(hEtaIntegral[det]!=0)delete hEtaIntegral[det];
//
//	if(tree!=0) delete tree;
}

bool TADCEventReader::SetTree(string fileName){//TTree *tree){
  pEvent=NULL;
	if(tree!=NULL) {tree->Delete();tree=NULL;}
	if(file!=NULL) {file->Delete();file=NULL;}
	tree=NULL;
	file=NULL;
//	cout<<"TADCEventReader-PATH: "<<sys->pwd()<<endl;
	std::cout<<"load File: \""<<fileName<<"\""<<endl;
//	stringstream fileString;
//	fileString<<sys->pwd()<<"/"<<fileName;
//	std::cout<<"Open "<<fileString.str()<<endl;
	file = (TFile*)TFile::Open(fileName.c_str());
	file->GetObject("tree",tree);
	if(tree==NULL)
	{
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
//	cout<<"TADCEventReader::isOK \""<<tree<<"\""<<!tree->IsZombie()<<endl;
	return (tree!=NULL&&!tree->IsZombie());
}

void TADCEventReader::SetBranchAddresses(){
  pEvent=0;
	//Event Header Branches
	if(tree->FindBranch("RunNumber")){
		tree->SetBranchAddress("RunNumber",&run_number);
		cout<<"Set Branch \"RunNumber\""<<endl;
	}
	else if(tree->FindBranch("runNumber")){
			tree->SetBranchAddress("runNumber",&run_number);
			cout<<"Set Branch \"runNumber\""<<endl;
	}
	if(tree->FindBranch("EventNumber")){
		tree->SetBranchAddress("EventNumber",&event_number);
		cout<<"Set Branch \"EventNumber\""<<endl;
	}
	if(tree->FindBranch("StoreThreshold")){
		tree->SetBranchAddress("StoreThreshold",&store_threshold);
		cout<<"Set Branch \"StoreThreshold\""<<endl;
	}
//	tree->SetBranchAddress("CMNEvent_flag",&CMNEvent_flag);

	if(tree->FindBranch("ZeroDivisorEvent_flag")){
		tree->SetBranchAddress("ZeroDivisorEvent_flag",&ZeroDivisorEvent_flag);
		cout<<"Set Branch \"ZeroDivisorEvent_flag\""<<endl;
	}//why do we have that????
	//Telescope Data Branches
	if(tree->FindBranch("D0X_NChannels")){
		tree->SetBranchAddress("D0X_NChannels",&Det_NChannels[0]);
		if(verbosity)cout<<"Set Branch \"D0X_NChannels\""<<endl;
	}
	if(tree->FindBranch("D0Y_NChannels")){
		tree->SetBranchAddress("D0Y_NChannels",&Det_NChannels[1]);
		if(verbosity)cout<<"Set Branch \"D0Y_NChannels\""<<endl;
		}
	if(tree->FindBranch("D1X_NChannels")){
		tree->SetBranchAddress("D1X_NChannels",&Det_NChannels[2]);
		if(verbosity)cout<<"Set Branch \"D1X_NChannels\""<<endl;
		}
	if(tree->FindBranch("D1Y_NChannels")){
		tree->SetBranchAddress("D1Y_NChannels",&Det_NChannels[3]);
		if(verbosity)cout<<"Set Branch \"D1Y_NChannels\""<<endl;
		}
	if(tree->FindBranch("D2X_NChannels")){
		tree->SetBranchAddress("D2X_NChannels",&Det_NChannels[4]);
		if(verbosity)cout<<"Set Branch \"D2X_NChannels\""<<endl;
		}
	if(tree->FindBranch("D2Y_NChannels")){
		tree->SetBranchAddress("D2Y_NChannels",&Det_NChannels[5]);
		if(verbosity)cout<<"Set Branch \"D2Y_NChannels\""<<endl;
		}
	if(tree->FindBranch("D3X_NChannels")){
		tree->SetBranchAddress("D3X_NChannels",&Det_NChannels[6]);
		if(verbosity)cout<<"Set Branch \"D3X_NChannels\""<<endl;
		}
	if(tree->FindBranch("D3Y_NChannels")){
		tree->SetBranchAddress("D3Y_NChannels",&Det_NChannels[7]);
		if(verbosity)cout<<"Set Branch \"D3Y_NChannels\""<<endl;
		}
	if(tree->FindBranch("Dia_NChannels")){
		tree->SetBranchAddress("Dia_NChannels",&Det_NChannels[8]);
		if(verbosity)cout<<"Set Branch \"Dia_NChannels\""<<endl;
		}
	if(tree->FindBranch("D0X_Channels")){
		tree->SetBranchAddress("D0X_Channels",&Det_Channels[0]);
		if(verbosity)cout<<"Set Branch \"D0X_Channels\""<<endl;
		}
	if(tree->FindBranch("D0Y_Channels")){
		tree->SetBranchAddress("D0Y_Channels",&Det_Channels[1]);
		if(verbosity)cout<<"Set Branch \"D0Y_Channels\""<<endl;
		}
	if(tree->FindBranch("D1X_Channels")){
		tree->SetBranchAddress("D1X_Channels",&Det_Channels[2]);
		if(verbosity)cout<<"Set Branch \"D1X_Channels\""<<endl;
		}
	if(tree->FindBranch("D1Y_Channels")){
		tree->SetBranchAddress("D1Y_Channels",&Det_Channels[3]);
		if(verbosity)cout<<"Set Branch \"D1Y_Channels\""<<endl;
		}
	if(tree->FindBranch("D2X_Channels")){
		tree->SetBranchAddress("D2X_Channels",&Det_Channels[4]);
		if(verbosity)cout<<"Set Branch \"D2X_Channels\""<<endl;
		}
	if(tree->FindBranch("D2Y_Channels")){
		tree->SetBranchAddress("D2Y_Channels",&Det_Channels[5]);
		if(verbosity)cout<<"Set Branch \"D2Y_Channels\""<<endl;
		}
	if(tree->FindBranch("D3X_Channels")){
		tree->SetBranchAddress("D3X_Channels",&Det_Channels[6]);
		if(verbosity)cout<<"Set Branch \"D3X_Channels\""<<endl;
		}
	if(tree->FindBranch("D3Y_Channels")){
		tree->SetBranchAddress("D3Y_Channels",&Det_Channels[7]);
		if(verbosity)cout<<"Set Branch \"D3Y_Channels\""<<endl;
		}
	if(tree->FindBranch("Dia_Channels")){
		tree->SetBranchAddress("Dia_Channels",&Det_Channels[8]);
		if(verbosity)cout<<"Set Branch \"Dia_Channels\""<<endl;
		}
	//tree->SetBranchAddress("Det_ADC",&Det_ADC[0][0]);
	if(tree->FindBranch("D0X_ADC")){
		tree->SetBranchAddress("D0X_ADC",&Det_ADC[0]);
		if(verbosity)cout<<"Set Branch \"D0X_ADC\""<<endl;
		}
	if(tree->FindBranch("D0Y_ADC")){
		tree->SetBranchAddress("D0Y_ADC",&Det_ADC[1]);
		if(verbosity)cout<<"Set Branch \"D0Y_ADC\""<<endl;
		}
	if(tree->FindBranch("D1X_ADC")){
		tree->SetBranchAddress("D1X_ADC",&Det_ADC[2]);
		if(verbosity)cout<<"Set Branch \"D1X_ADC\""<<endl;
		}
	if(tree->FindBranch("D1Y_ADC")){
		tree->SetBranchAddress("D1Y_ADC",&Det_ADC[3]);
		if(verbosity)cout<<"Set Branch \"D1Y_ADC\""<<endl;
		}
	if(tree->FindBranch("D2X_ADC")){
		tree->SetBranchAddress("D2X_ADC",&Det_ADC[4]);
		if(verbosity)cout<<"Set Branch \"D2X_ADC\""<<endl;
		}
	if(tree->FindBranch("D2Y_ADC")){
		tree->SetBranchAddress("D2Y_ADC",&Det_ADC[5]);
		if(verbosity)cout<<"Set Branch \"D2Y_ADC\""<<endl;
		}
	if(tree->FindBranch("D3X_ADC")){
		tree->SetBranchAddress("D3X_ADC",&Det_ADC[6]);
		if(verbosity)cout<<"Set Branch \"D3X_ADC\""<<endl;
		}
	if(tree->FindBranch("D3Y_ADC")){
		tree->SetBranchAddress("D3Y_ADC",&Det_ADC[7]);
		if(verbosity)cout<<"Set Branch \"D3Y_ADC\""<<endl;
		}
	//tree->SetBranchAddress("Dia_ADC",&Dia_ADC);
	if(tree->FindBranch("DiaADC")){
		tree->SetBranchAddress("DiaADC",&Dia_ADC);
		if(verbosity)cout<<"Set Branch \"DiaADC\""<<endl;
		}
	if(tree->FindBranch("D0X_PedMean")){
		tree->SetBranchAddress("D0X_PedMean",&Det_PedMean[0]);
		if(verbosity)cout<<"Set Branch \"D0X_PedMean\""<<endl;
		}
	if(tree->FindBranch("D0Y_PedMean")){
		tree->SetBranchAddress("D0Y_PedMean",&Det_PedMean[1]);
		if(verbosity)cout<<"Set Branch \"D0Y_PedMean\""<<endl;
		}
	if(tree->FindBranch("D1X_PedMean")){
		tree->SetBranchAddress("D1X_PedMean",&Det_PedMean[2]);
		if(verbosity)cout<<"Set Branch \"D1X_PedMean\""<<endl;
	}
	if(tree->FindBranch("D1Y_PedMean")){
		tree->SetBranchAddress("D1Y_PedMean",&Det_PedMean[3]);
		if(verbosity)cout<<"Set Branch \"D1Y_PedMean\""<<endl;
	}
	if(tree->FindBranch("D2X_PedMean")){
		tree->SetBranchAddress("D2X_PedMean",&Det_PedMean[4]);
		if(verbosity)cout<<"Set Branch \"D2X_PedMean\""<<endl;
	}
	if(tree->FindBranch("D2Y_PedMean")){
		tree->SetBranchAddress("Dia_Channels",&Det_PedMean[5]);
		if(verbosity)cout<<"Set Branch \"Dia_Channels\""<<endl;
		}
	if(tree->FindBranch("D3X_PedMean")){
		tree->SetBranchAddress("D3X_PedMean",&Det_PedMean[6]);
		if(verbosity)cout<<"Set Branch \"D3X_PedMean\""<<endl;
		}
	if(tree->FindBranch("D3Y_PedMean")){
		tree->SetBranchAddress("D3Y_PedMean",&Det_PedMean[7]);
		if(verbosity)cout<<"Set Branch \"D3Y_PedMean\""<<endl;
		}
	if(tree->FindBranch("Dia_PedMean")){
		tree->SetBranchAddress("Dia_PedMean",&Det_PedMean[8]);
		if(verbosity)cout<<"Set Branch \"Dia_PedMean\""<<endl;
		}
	if(tree->FindBranch("D0X_PedWidth")){
		tree->SetBranchAddress("D0X_PedWidth",&Det_PedWidth[0]);
		if(verbosity)cout<<"Set Branch \"D0X_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D0Y_PedWidth")){
		tree->SetBranchAddress("D0Y_PedWidth",&Det_PedWidth[1]);
		if(verbosity)cout<<"Set Branch \"D0Y_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D1X_PedWidth")){
		tree->SetBranchAddress("D1X_PedWidth",&Det_PedWidth[2]);
		if(verbosity)cout<<"Set Branch \"D1X_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D1Y_PedWidth")){
		tree->SetBranchAddress("D1Y_PedWidth",&Det_PedWidth[3]);
		if(verbosity)cout<<"Set Branch \"D1Y_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D2X_PedWidth")){
		tree->SetBranchAddress("D2X_PedWidth",&Det_PedWidth[4]);
		if(verbosity)cout<<"Set Branch \"D2X_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D2Y_PedWidth")){
		tree->SetBranchAddress("D2Y_PedWidth",&Det_PedWidth[5]);
		if(verbosity)cout<<"Set Branch \"D2Y_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D3X_PedWidth")){
		tree->SetBranchAddress("D3X_PedWidth",&Det_PedWidth[6]);
		if(verbosity)cout<<"Set Branch \"D3X_PedWidth\""<<endl;
		}
	if(tree->FindBranch("D3Y_PedWidth")){
		tree->SetBranchAddress("D3Y_PedWidth",&Det_PedWidth[7]);
		if(verbosity)cout<<"Set Branch \"D3Y_PedWidth\""<<endl;
		}
	if(tree->FindBranch("Dia_PedWidth")){
		tree->SetBranchAddress("Dia_PedWidth",&Det_PedWidth[8]);
		if(verbosity)cout<<"Set Branch \"Dia_PedWidth\""<<endl;
	}
	if(tree->FindBranch("PedestalMean")){
		tree->SetBranchAddress("PedestalMean",&pedestalMean);
		if(verbosity)cout<<"Set Branch \"PedestalMean\""<<endl;
		}
	if(tree->FindBranch("PedestalSigma")){
		tree->SetBranchAddress("PedestalSigma",&pedestalSigma);
		if(verbosity)cout<<"Set Branch \"PedestalSigma\""<<endl;
		}
	if(tree->FindBranch("diaPedestalMean")){
	    tree->SetBranchAddress("diaPedestalMean",&diaPedestalMean);
	    if(verbosity)cout<<"Set Branch \"diaPedestalMean\""<<endl;
	    }
	  if(tree->FindBranch("diaPedestalSigma")){
	    tree->SetBranchAddress("diaPedestalSigma",&diaPedestalSigma);
	    if(verbosity)cout<<"Set Branch \"diaPedestalSigma\""<<endl;
	    }
	  if(tree->FindBranch("diaPedestalMeanCMN")){
	    tree->SetBranchAddress("diaPedestalMeanCMN",&diaPedestalMeanCMN);
	    if(verbosity)cout<<"Set Branch \"diaPedestalMeanCMN\""<<endl;
	  }
	  if(tree->FindBranch("diaPedestalSigmaCMN")){
	    tree->SetBranchAddress("diaPedestalSigmaCMN",&diaPedestalSigmaCMN);
	    if(verbosity)cout<<"Set Branch \"diaPedestalSigmaCMN\""<<endl;
	  }
//	if(tree->FindBranch("clusters")){
//		//tree->SetBranchAddress("clusters",&pVecvecCluster);
//		if(verbosity)cout<<"Set Branch \"clusters\""<<endl;
//		}
	if(tree->FindBranch("isDetMasked")){
		tree->SetBranchAddress(	"isDetMasked",&bIsDetMasked);
		if(verbosity)cout<<"Set Branch \"isDetMasked\""<<endl;
		}
	if(tree->FindBranch("hasValidSiliconTrack")){
		tree->SetBranchAddress("hasValidSiliconTrack",&bValidSiliconTrack);
		if(verbosity)cout<<"Set Branch \"hasValidSiliconTrack\""<<endl;
		}
	if(tree->FindBranch("nDiamondHits")){
		tree->SetBranchAddress("nDiamondHits",&nDiamondClusters);
		if(verbosity)cout<<"Set Branch \"nDiamondHits\""<<endl;
		}
	if(tree->FindBranch("isInFiducialCut")){
		tree->SetBranchAddress("isInFiducialCut",&bIsInFiducialCut);
		if(verbosity)cout<<"Set Branch \"isInFiducialCut\""<<endl;
	}
//	if(tree->FindBranch("isDiaMasked")){
//		tree->SetBranchAddress("isDiaMasked",&this->maskedDiaClusters);
//		if(verbosity)cout<<"Set Branch \"isDiaMasked\""<<endl;
//	}
	if(tree->FindBranch("event")){
		tree->SetBranchAddress("event",&pEvent);
		if(verbosity)cout<<"Set Branch \"event\""<<endl;
	}
	else
		if(verbosity)cout<<" \"event\" not found..."<<endl;
	if(tree->FindBranch("useForAlignment")){
		tree->SetBranchAddress("useForAlignment",&this->bUseForAlignment);
		if(verbosity)cout<<"Set Branch \"useForAlignment\""<<endl;
	}
	if(tree->FindBranch("useForAnalysis")){
		tree->SetBranchAddress("useForAnalysis",&this->bUseForAnalysis);
		if(verbosity)cout<<"Set Branch \"useForAnalysis\""<<endl;
	}
	if(tree->FindBranch("useForAlignment")){
		tree->SetBranchAddress("useForSiliconAlignment",&this->bUseForAlignment);
		if(verbosity)cout<<"Set Branch \"useForSiliconAlignment\""<<endl;
	}

	if(tree->FindBranch("cmnCorrection")){
	  tree->SetBranchAddress("cmnCorrection",&this->bCMNoiseCorrected);
	  cout<<"Set Branch \n isCMNCOrrected\""<<endl;
	}
	if(tree->FindBranch("commonModeNoise")){
	  tree->SetBranchAddress("commonModeNoise",&this->cmNoise);
	  cout<<"SET BRANCH ADDRESS: \"CMN Noise\""<<endl;
	}

//	vector<bool> isDiaMasked;//thediamond plane contains a cluster wit a masked channel (size of nDiamondHits)
//	UInt_t nDiamondHits; //number of clusters in diamond plane;
	cout<<"DONE"<<endl;

}

void TADCEventReader::initialiseTree(){
  bCMNoiseCorrected=false;
	cout<<"initialise tree with "<<tree->GetEntries()<<" Entires."<<endl;
	current_event = 1;
	tree->GetEvent(current_event);
	cout<< "Loaded first event in Tree: "<<event_number<<endl;
	cout<< "RunNumber is: "<<run_number<<endl;
	cout<< "StoreThreshold is: "<<store_threshold<<endl;
}

bool TADCEventReader::GetNextEvent(){
	return LoadEvent(current_event+1);
}
bool TADCEventReader::LoadEvent(UInt_t EventNumber){
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
	if(det<TPlaneProperties::getNSiliconDetectors()&& TPlaneProperties::getNChannels(det)>ch)
		return Det_ADC[det][ch];
	return -1;
}

UChar_t TADCEventReader::getDet_Channels(UInt_t i , UInt_t j) const
{
	if (i>8||j>255){
		cout<<"TADCEventReader::getDet_Channels not Valid "<<i<<" "<<j<<endl;
		exit (-1);
	}
    return Det_Channels[i][j];
}

UInt_t TADCEventReader::getDet_NChannels(UInt_t det)  const
{
    return Det_NChannels[det];
}

Float_t TADCEventReader::getDet_PedMean(UInt_t i, UInt_t j) const
{
    return Det_PedMean[i][j];
}

Float_t TADCEventReader::getDet_PedWidth(UInt_t i, UInt_t j) const
{
    return Det_PedWidth[i][j];
}

Int_t TADCEventReader::getDia_ADC(UInt_t ch) {
  if(ch<TPlaneProperties::getNChannelsDiamond())
    return (Int_t)Dia_ADC[ch];
  return -1;
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


TFile *TADCEventReader::getFile() const
{
	return this->file;
}

/**
 * give the pedestal value of det and ch
 * if cmnCorrected and det == diamond detector
 *   it will add CMN and pedestalValue
 * @returns pedestal or pedestal+CMN
 */
Float_t TADCEventReader::getPedestalMean(UInt_t det, UInt_t ch,bool cmnCorrected)
{
	if(TPlaneProperties::isSiliconDetector(det)&&ch<TPlaneProperties::getNChannels(det))
		return this->pedestalMean[det][ch];
	if(TPlaneProperties::isDiamondDetector(det))
	  return getDiaPedestalMean(ch,cmnCorrected);
	return -99999;
}

Float_t TADCEventReader::getPedestalSigma(UInt_t det, UInt_t ch,bool cmnCorrected)
{
	if(TPlaneProperties::isSiliconDetector(det)&&ch<TPlaneProperties::getNChannels(det))
		if(this->pedestalSigma[det][ch]>=0)
			return this->pedestalSigma[det][ch];
	if(TPlaneProperties::isDiamondDetector(det))
	  return getDiaPedestalSigma(ch,cmnCorrected);
	return 0;
}


Float_t TADCEventReader::getDiaPedestalMean(UInt_t ch, bool cmnCorrected){
  if(ch<TPlaneProperties::getNChannelsDiamond()){
    if(cmnCorrected)
      return this->diaPedestalMeanCMN[ch];
    else
      return this->diaPedestalMean[ch];
  }
  return -99999;
}

Float_t TADCEventReader::getDiaPedestalSigma(UInt_t ch, bool cmnCorrected){
  if(ch<TPlaneProperties::getNChannelsDiamond()){
    if(cmnCorrected)
      return this->diaPedestalSigmaCMN[ch];
    else
      return this->diaPedestalSigma[ch];
  }
  return 0;
}

Int_t TADCEventReader::getAdcValue(UInt_t det,UInt_t ch){
	if(TPlaneProperties::isSiliconDetector(det))
		return (Int_t)this->getDet_ADC(det,ch);
	if (TPlaneProperties::isDiamondDetector(det))
	    return (Int_t)this->getDia_ADC(ch);
	return -1;
}

Float_t TADCEventReader::getSignalInSigma(UInt_t det, UInt_t ch, bool cmnCorrected)
{
	if(getPedestalSigma(det,ch,cmnCorrected)==0)
		return 0;
	else
		return (this->getSignal(det,ch,cmnCorrected)/this->getPedestalSigma(det,ch,cmnCorrected));
}

TCluster TADCEventReader::getCluster(UInt_t det, UInt_t cl)
{
	if(pEvent!=NULL)
		return this->pEvent->getCluster(det,cl);
	return TCluster();
}
TCluster TADCEventReader::getCluster(UInt_t plane,TPlaneProperties::enumCoordinate cor, UInt_t cl){
	if(pEvent!=NULL)
	return this->pEvent->getCluster(plane,cor,cl);
	return TCluster();
}

UInt_t TADCEventReader::getClusterSize(UInt_t det,UInt_t cl)
{
	if(pEvent!=NULL)
	return pEvent->getClusterSize(det,cl);
	return 0;
}
UInt_t TADCEventReader::getClusterSeedSize(UInt_t det,UInt_t cl)
{
	if(pEvent!=NULL)
	return pEvent->getClusterSeedSize(det,cl);
	return 0;
}

void TADCEventReader::checkADC(){
	this->LoadEvent(100);
	for(int ch=0;ch<256;ch++)
		cout<<this->getAdcValue(0,ch)<<" "<<this->getAdcValue(1,ch)<<" "<<this->getAdcValue(8,ch)<<endl;
}

bool TADCEventReader::isSaturated(UInt_t det, UInt_t ch)
{
	if(det<TPlaneProperties::getNDetectors())
		return getAdcValue(det,ch)>=TPlaneProperties::getMaxSignalHeight(det);
	return true;
}

Float_t TADCEventReader::getRawSignal(UInt_t det, UInt_t ch,bool cmnCorrected){
  if(det>=9)return -9999999;
  Float_t cmn = getCMNoise();
  Int_t adc = getAdcValue(det,ch);
//  float ped = getPedestalMean(det,ch,false);
//  float pedCMN = getPedestalMean(det,ch,true);
  Float_t pedReal= getPedestalMean(det,ch,cmnCorrected);
  if (!cmnCorrected||TPlaneProperties::isSiliconDetector(det))
         cmn=0;
//  if (TPlaneProperties::isDiamondDetector(det)){
////    cout<<getEvent_number(  )<<" "<<det<<" "<<ch<<" "<<cmnCorrected<<endl;
////    cout<<"raw signal det:"<<det<<" ch:"<<ch<<": "<<adc<<"-"<<ped<<"="<<adc-ped<<endl;
////    cout<<"               "<<det<<" ch:"<<ch<<": "<<adc<<"-"<<pedCMN<<"-"<<cmn<<"="<<adc-pedCMN-cmn<<endl;
//
//    retVal = adc-pedReal-cmn;
////    cout<<"               "<<det<<" ch:"<<ch<<": "<<adc<<"-"<<pedReal<<"+"<<cmn<<"="<<retVal<<endl;
//  }
//  else
//    cmn=0;
 Float_t retVal = adc-pedReal-cmn;
  return retVal;
}

Float_t TADCEventReader::getRawSignalInSigma(UInt_t det, UInt_t ch,bool cmnCorrected){
  if(det>=TPlaneProperties::getNDetectors()||(getPedestalSigma(det,ch,cmnCorrected)<=0))
	  return -99999999;
  return (getRawSignal(det,ch,cmnCorrected)/getPedestalSigma(det,ch,cmnCorrected));
}

Float_t TADCEventReader::getSignal(UInt_t det, UInt_t ch,bool cmnCorrected)
{
	if(det>=9) return -9999999;

	Float_t signal =getRawSignal(det,ch,cmnCorrected);
	if(signal<0)return 0;

	return signal;
}

UInt_t TADCEventReader::getNClusters(UInt_t det)
{
	if(det<TPlaneProperties::getNDetectors()){

		UInt_t nClusters = this->pEvent->getNClusters(det);
		if(verbosity>1){
			cout<<"TADCEventReader::getNClusters of det "<<det<<": "<<nClusters<<endl;
			pEvent->setVerbosity(verbosity);
		}
		return nClusters;
	}
	return 0;
}

bool TADCEventReader::isValidTrack()
{
	return this->bValidSiliconTrack; // one & only one hit in silicone planes
}

UInt_t TADCEventReader::getNDiamondClusters()
{
	return this->nDiamondClusters;
}

bool TADCEventReader::isInFiducialCut()
{
	return this->bIsInFiducialCut;
}


bool TADCEventReader::isDetMasked()
{
	return this->bIsDetMasked;
}

TEvent* TADCEventReader::getEvent()
{
	return this->pEvent;
}

void TADCEventReader::setVerbosity(UInt_t verbosity)
{
	this->verbosity=verbosity;
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
			cout<<"name of it: "<<obj->GetName()<<endl;
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

TH1F *TADCEventReader::getEtaIntegral(UInt_t det)
{
	if(!bEtaIntegrals||det>8)
		return 0;
	return hEtaIntegral[det];
}

void TADCEventReader::LoadEtaDistributions(UInt_t runNumber){
	bEtaIntegrals=true;
	stringstream etaFileName;

  if(!TSettings::existsDirectory(etaDistributionPath))
    etaFileName<<"etaCorrection."<<runNumber<<".root";
  else
    etaFileName<<etaDistributionPath;

	TFile *fEtaDis = TFile::Open(etaFileName.str().c_str());
	if(fEtaDis==0){
		cout<<"EtaDistribution File \""<<etaFileName.str()<<"\" do not exist"<<endl;
		bEtaIntegrals=false;
		if(etaDistributionPath.size()!=0){char t; cin>>t;}
		return;
	}
	cout<<etaFileName.str()<<endl;
//	char t; cin>>t;
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		stringstream objectName;
		objectName<<"hEtaIntegral_"<<det;
		TH1F *histo = 0;
		fEtaDis->GetObject(objectName.str().c_str(),histo);
		file->cd();
		bEtaIntegrals=bEtaIntegrals&&(histo!=0);
		if(histo)
			hEtaIntegral[det]=(TH1F*)histo->Clone();
		else
			cerr<<"Object \""<<objectName.str()<<"\" does not exist"<<endl;
	}

}

void TADCEventReader::setEtaDistributionPath(std::string path)
{
  etaDistributionPath=path;
}


