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
  initialiseResults();
}

TResults::TResults(TSettings *settings) {
  path = gSystem->pwd();
  initialiseResults();
  cout<<"New TResults with settings "<<settings<<"\t run: "<<settings->getRunNumber()<<endl;
  openResults(settings);
//  this->settings=settings;
}

TResults::~TResults() {
  // TODO Auto-generated destructor stub
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
}

void TResults::initialiseResults(){
  seedSigma.resize(TPlaneProperties::getNDetectors(),-9999);
  hitSigma.resize(TPlaneProperties::getNDetectors(),-9999);
  noise.resize(TPlaneProperties::getNDetectors(),-9999);
}


void TResults::openResults(TSettings *settings){
  cout<<"open Results with settings "<<settings<<"\t run: "<<settings->getRunNumber()<<flush;
//  this->Settings = *settings;
  runnumber = settings->getRunNumber();
  cout<<" "<<runnumber<<endl;
  std::stringstream resultsFile;
  resultsFile<<path<<"/"<<runnumber<<"/Results."<<runnumber<<".root";
  cout<<resultsFile.str()<<endl;
  TFile *file =  new TFile(resultsFile.str().c_str(),"READ");
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
    oldResults = (TResults*)file->Get(name.str().c_str());
    cout<<"old Results: "<<oldResults<<endl;

    if(oldResults==0){
      cerr<< "Something is wrong, results does not exists..."<<endl;
      return;
    }
    cout<<oldResults->IsZombie()<<endl;
//    cout<<oldResults->IsA()->ClassName()<<endl;
    oldResults->getLastUpdateDate().Print();
    cout<<"LAST UPDATE ON "<<oldResults->getLastUpdateDate().AsString()<<endl;
    this->inheritOldResults(*oldResults);
    for(UInt_t det=0;det<TPlaneProperties::getNDetectors();det++){
      seedSigma.at(det)=settings->getClusterSeedFactor(det);
      hitSigma.at(det)=settings->getClusterHitFactor(det);
    }

  }
  Print();

}

void TResults::saveResults(){
  lastUpdate=TDatime();
  std::stringstream fileName;
  fileName<<path<<"/"<<runnumber<<"/Results."<<runnumber<<".root";
  std::stringstream name;
  name <<"results_"<<runnumber;
  this->SetName(name.str().c_str());
  TFile *file =  new TFile(fileName.str().c_str(),"RECREATE");
  file->cd();
  this->Write();
  file->Close();
}

void TResults::Print(){
  getLastUpdateDate().Print();
  cout<<getLastUpdateDate().AsString()<<endl;
  cout<<"det\tseed hit  noise"<<endl;
  for(UInt_t det=0; det<TPlaneProperties::getNDetectors();det++){
    cout<<det<<"\t"<<std::setw(4)<<this->seedSigma[det]<<" "<<std::setw(4)<<this->hitSigma[det]<<" "<<std::setw(4)<<this->noise[det]<<endl;
  }
//  settings->Print();
}

void TResults::SetNoise(UInt_t det,Float_t detNoise){
  if(det>=TPlaneProperties::getNDetectors())
    return;
  if(noise.size()<=det)
    noise.resize(det+1,0);
  noise[det]=detNoise;
  cout<<"Set Results: Noise of det "<<det<<": "<<detNoise<<endl;
}



void TResults::setAlignment(TDetectorAlignment newAlignment)
{
  alignment= newAlignment;
}
