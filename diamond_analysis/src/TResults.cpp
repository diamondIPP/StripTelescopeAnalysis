/*
 * TResults.cpp
 *
 *  Created on: May 29, 2012
 *      Author: bachmair
 */

#include "../include/TResults.hh"

ClassImp(TResults);

using namespace std;

TResults::TResults(UInt_t runnumber){
	path = gSystem->pwd();
	this->runnumber = runnumber;
	TString name =  TString::Format("results_%d",this->runnumber);
	this->SetName(name);
	initialiseResults();
	rootFileName = TString::Format("results.%d.root",this->runnumber);
}

TResults::TResults(TSettings *settings) {
	path = gSystem->pwd();
	initialiseResults();
	cout<<"New TResults with settings "<<settings;
	if (settings) cout << "\t run: "<<settings->getRunNumber()<<endl;
	else cout<<"\t run: XXXXXXXX "<<endl;
	openResults(settings);
	textFileName = settings->getAbsoluteOuputPath(true);
	textFileName.Append(TString::Format("/results_%d.txt",this->runnumber));
	if (settings)
		rootFileName = settings->getResultsRootFilePath();
	else
		rootFileName = TString::Format("results.%d.root",this->runnumber);
	//  this->settings=settings;
}

TResults::~TResults() {
	// TODO Auto-generated destructor stub
	cout<<"delete Results of run "<<runnumber<<endl;
}

TResults::TResults(const   TResults& rhs){//copy constructor
	initialiseResults();
	inheritOldResults(rhs);
}
//
//TResults::TResults &operator=(const   TResults &src){ //class assignment function
//
//}

void TResults::inheritOldResults(const TResults & rhs)
{

	this->seedSigma.clear();
	for(UInt_t det=0;det<rhs.seedSigma.size();det++)this->seedSigma.push_back(rhs.seedSigma[det]);
	this->hitSigma.clear();
	for(UInt_t det=0;det<rhs.hitSigma.size();det++)this->hitSigma.push_back(rhs.hitSigma[det]);
	this->noise.clear();
	for(UInt_t det=0;det<rhs.noise.size();det++)this->noise.push_back(rhs.noise[det]);
	diaCMCNoise = rhs.diaCMCNoise;
	CMN = rhs.CMN;
	mean2outOf10_normal = rhs.mean2outOf10_normal;
	mean2outOf10_trans = rhs.mean2outOf10_trans;
	meanNoutOfN_normal.clear();
	for(UInt_t cl=0;cl<rhs.meanNoutOfN_normal.size();cl++)this->meanNoutOfN_normal.push_back(rhs.meanNoutOfN_normal[cl]);
	meanNoutOfN_trans.clear();
	for(UInt_t cl=0;cl<rhs.meanNoutOfN_trans.size();cl++)this->meanNoutOfN_trans.push_back(rhs.meanNoutOfN_trans[cl]);

	mp2outOf10_normal = rhs.mp2outOf10_normal;
	mp2outOf10_trans = rhs.mp2outOf10_trans;
	width2outOf10_normal = rhs.width2outOf10_normal;
	width2outOf10_trans = rhs.width2outOf10_trans;
	gSigma2outOf10_normal = rhs.gSigma2outOf10_normal;
	gSigma2outOf10_trans = rhs.gSigma2outOf10_trans;
	this->lastUpdate = rhs.lastUpdate;
}

void TResults::initialiseResults(){
	seedSigma.resize(TPlaneProperties::getNDetectors(),-1);
	hitSigma.resize(TPlaneProperties::getNDetectors(),-1);
	noise.resize(TPlaneProperties::getNDetectors(),-1);
	meanNoutOfN_normal.resize(TPlaneProperties::getMaxTransparentClusterSize(TPlaneProperties::getDetDiamond()),-1);
	meanNoutOfN_trans.resize(TPlaneProperties::getMaxTransparentClusterSize(TPlaneProperties::getDetDiamond()),-1);
	doubleGaus1_normal = -1;
	doubleGaus2_normal = -1;
	doubleGaus1_trans = -1;
	doubleGaus2_trans = -1;

	singleGausShort_normal = -1;
	singleGausShort_trans = -1;
	singleGaus_normal = -1;
	singleGaus_trans = -1;
}


void TResults::openResults(TSettings *settings){
	cout<< settings<<endl;
	//	if (!settings)return;
	cout<<"open Results with settings "<<settings<<"\t run: "<<settings->getRunNumber()<<flush;
	rootFileName = settings->getResultsRootFilePath();
	//  this->Settings = *settings;
	runnumber = settings->getRunNumber();
	cout<<" "<<runnumber<<endl;
	//	std::stringstream resultsFile;
	//	resultsFile<<path<<"/"<<runnumber<<"/Results."<<runnumber<<".root";
	cout<<((string)(rootFileName))<<endl;

	TFile *file =  new TFile(rootFileName,"READ");
	TResults *oldResults;

	if(file->IsZombie()) {
		cout << "FIle does not exists, create new File!"<<endl;
		delete file;
		oldResults = new TResults(runnumber);
	}
	else{
		file->GetListOfKeys()->Print();
		stringstream name;
		name << "results_"<<runnumber;
		cout<<"Name of key: \""<<name.str()<<"\""<<endl;
		TIter next(file->GetListOfKeys());
		TKey *key;
		oldResults = 0;
		bool foundResult = false;
		while ((key=(TKey*)next())) {
			TString className = key->GetClassName();
			if (className.Contains("TResults")){
				if (foundResult = true && oldResults){
					cout <<" found a second Results key in the file..."<<oldResults->GetName()<<" "<<key->GetName()<<endl;
				}
				else{
					cout<<"FOUND a valid class: "<<key->GetName()<<endl;
					oldResults =  (TResults*)key->ReadObj();
					foundResult = true;
				}
			}
		}
//		oldResults = (TResults*)file->FindGet(name.str().c_str());
		cout<<"old Results: "<<oldResults<<flush;
		if (oldResults) cout<<" "<<oldResults->GetName()<<flush;

		if(oldResults==0){
			cerr<< "Something is wrong, results does not exists..."<<endl;
			return;
		}
		cout<<oldResults->IsZombie()<<endl;
//		cout<<"LAST UPDATE ON "<<oldResults->getLastUpdateDate().AsString()<<endl;
		this->inheritOldResults(*oldResults);
		setResultsFromSettings(settings);

	}
	Print();
}

void TResults::setResultsFromSettings(TSettings* settings){
	for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
		seedSigma.at(det)=settings->getClusterSeedFactor(det,0);
		hitSigma.at(det)=settings->getClusterHitFactor(det,0);
	}
}

void TResults::saveResults(TString name){
	lastUpdate=TDatime();
	//	std::stringstream fileName;
	//	fileName<<path<<"/"<<runnumber<<"/Results."<<runnumber<<".root";
	//	cout<<"save File:"<<fileName.str()<<endl;
	//	std::stringstream name;
	//	name <<"results_"<<runnumber;
	//	this->SetName(name.str().c_str());
	//	TFile *file =  new TFile(fileName.str().c_str(),"RECREATE");
	//	file->cd();
	//	this->Write();
	//	this->Write(fileName.str().c_str());
	cout<< "Save "<<name<<endl;
	//	this->Write(name);
	//	TFile *file =  new TFile(name.Append(".bak"),"RECREATE");
	//	this->Write();
	//	file->Close();
}

void TResults::Print(){
	//	getLastUpdateDate().Print();
	//	cout<<"\t"<<ÊgetLastUpdateDate().AsString()<<endl;
	cout<<getLastUpdateDate().AsString()<<endl;
	cout<<"\tdet\tseed  \thit  \tnoise"<<endl;
	for(UInt_t det=0; det<TPlaneProperties::getNDetectors();det++){
		cout<<"\t"<<det<<"\t"<<std::setw(4)<<this->seedSigma[det]<<"\t"<<std::setw(4)<<this->hitSigma[det]<<" \t"<<std::setw(4)<<this->noise[det]<<endl;
	}
	//	alignment.Print();

	//  settings->Print();z
	cout << "\tmean2outOf10_normal   = " << mean2outOf10_normal << "\n"
			<< "\tmp2outOf10_normal     = " << mp2outOf10_normal << "\n"
			<< "\tmp2outOf10_normal     = " << mp2outOf10_normal << "\n"
			<< "\tgSigma2outOf10_normal = " << gSigma2outOf10_normal<<endl;
	cout<<"DONE"<<endl;
	createOutputTextFile();
}

void TResults::SetNoise(UInt_t det,Float_t detNoise){
	if(det>=TPlaneProperties::getNDetectors())
		return;
	if(noise.size()<=det)
		noise.resize(det+1,0);
	noise[det]=detNoise;
	cout<<"Set Results: Noise of det "<<det<<": "<<detNoise<<endl;
}


void TResults::setAlignment(TDetectorAlignment* newAlignment)
{
	this->alignment= *newAlignment;
	cout<<"TResults:SetAlignment"<<endl;
}


int TResults::getRunNumber() const {
	return runnumber;
}

void TResults::setRunNumber(UInt_t runNumber) {
	this->runnumber = runNumber;
}


void TResults::setDoubleGaussianResolution(Float_t gaus1,Float_t gaus2, TSettings::alignmentMode mode){
	if (gaus2 > gaus1){
		Float_t val = gaus1;
		gaus1 = gaus2;
		gaus2 = val;
	}
	if (mode == TSettings::normalMode){
		doubleGaus1_normal = gaus1;
		doubleGaus2_normal = gaus2;
	}
	else if ( mode == TSettings::transparentMode){
		doubleGaus1_trans = gaus1;
		doubleGaus2_trans = gaus2;
	}
}

void TResults::setSingleGaussianResolution(Float_t gaus,TSettings::alignmentMode mode){
	if (mode == TSettings::normalMode){
		singleGaus_normal = gaus;
	}
	else if ( mode == TSettings::transparentMode){
		singleGaus_trans = gaus;
	}
}
void TResults::setSingleGaussianShortResolution(Float_t gaus,TSettings::alignmentMode mode){
	if (mode == TSettings::normalMode){
		singleGausShort_normal = gaus;
	}
	else if ( mode == TSettings::transparentMode){
		singleGausShort_trans = gaus;
	}
}

void TResults::setPH_2outOf10(Float_t mean, Float_t mp, Float_t width, Float_t gSigma, TSettings::alignmentMode mode){
	if (mode == TSettings::normalMode){
		mean2outOf10_normal = mean;
		mp2outOf10_normal = mp;
		width2outOf10_normal = width;
		gSigma2outOf10_normal = gSigma;
	}
	else if ( mode == TSettings::transparentMode){
		mean2outOf10_trans = mean;
		mp2outOf10_trans = mp;
		width2outOf10_trans = width;
		gSigma2outOf10_trans = gSigma;
	}
}

void TResults::setPH_NoutOfN(vector<Float_t> vecPHNoutOfN, TSettings::alignmentMode mode){
	if(meanNoutOfN_trans.size()<vecPHNoutOfN.size())
		meanNoutOfN_trans.resize(vecPHNoutOfN.size(),-1);
	if(meanNoutOfN_normal.size()<vecPHNoutOfN.size())
		meanNoutOfN_normal.resize(vecPHNoutOfN.size(),-1);

	if(mode == TSettings::normalMode){
		for (UInt_t i = 0; i< vecPHNoutOfN.size();i++)
			meanNoutOfN_normal.at(i) = vecPHNoutOfN.at(i);
	}
	else if (mode == TSettings::transparentMode){
		for (UInt_t i = 0; i< vecPHNoutOfN.size();i++)
			meanNoutOfN_trans.at(i) = vecPHNoutOfN.at(i);
	}

}

void TResults::setCMNoise(Float_t cmn){
	this->CMN = cmn;
}

void TResults::createOutputTextFile(){
//	cout<<"CREATE OUPTUT TEXT FILE"<<endl;
	ofstream myfile;
	myfile.open (textFileName, ios::out	|ios::trunc);
	UInt_t det  = 8;
	myfile << "#\t     \t      \tCMN\t";
	myfile << "normal\t      \t     \t      \t    \t    \t";
	myfile << "transAlign    \t     \t       \t    \t   \t";
	myfile << "Res. normal\t\t   \t";
	myfile << "Res. trans\t\t   \t";
	myfile << endl;
	myfile << "#RUN\tnoise\tCMCnoi\tCMN\t";
	myfile << "m2/10\tmp2/10\tw2/10\tsig2/10\tm2/2\tm4\4\t";
	myfile << "m2/10\tmp2/10\tw2/10\tsig2/10\tm2/2\tm4\4\t";
	myfile << "res1\tres2\tres3\tres4\t";
	myfile << "res1\tres2\tres3\tres4";
	myfile << endl;
	myfile << " ";
	myfile << TString::Format("%6d\t",runnumber);
	myfile << TString::Format("%2.2f\t",noise[det]);
	myfile << TString::Format("%2.2f\t",CMN);
	myfile << TString::Format("%2.2f\t",diaCMCNoise);
//	myfile <<" \t";
	myfile << TString::Format("%2.1f\t",mean2outOf10_normal);
	myfile << TString::Format("%2.1f\t",mp2outOf10_normal);
	myfile << TString::Format("%2.2f\t",width2outOf10_normal);
	myfile << TString::Format("%2.2f\t",gSigma2outOf10_normal);
	myfile << TString::Format("%2.1f\t",meanNoutOfN_normal[1]);
	myfile << TString::Format("%2.1f\t",meanNoutOfN_normal[3]);
//	myfile <<" \t";
	myfile << TString::Format("%2.1f\t",mean2outOf10_trans);
	myfile << TString::Format("%2.1f\t",mp2outOf10_trans);
	myfile << TString::Format("%2.2f\t",width2outOf10_trans);
	myfile << TString::Format("%2.2f\t",gSigma2outOf10_trans);
	myfile << TString::Format("%2.1f\t",meanNoutOfN_trans[1]);
	myfile << TString::Format("%2.1f\t",meanNoutOfN_trans[3]);
	//	myfile <<" \t";
	myfile << TString::Format("%2.2f\t",doubleGaus1_normal);
	myfile << TString::Format("%2.2f\t",doubleGaus2_normal);
	myfile << TString::Format("%2.2f\t",singleGausShort_normal);
	myfile << TString::Format("%2.2f\t",singleGaus_normal);
	//	myfile <<" \t";
	myfile << TString::Format("%2.2f\t",doubleGaus1_trans);
	myfile << TString::Format("%2.2f\t",doubleGaus2_normal);
	myfile << TString::Format("%2.2f\t",singleGausShort_trans);
	myfile << TString::Format("%2.2f\t",singleGaus_trans);
	myfile << endl;

	myfile.close();
//	cout<<"CREATE OUPTUT TEXT FILE DONE "<<endl;
}
