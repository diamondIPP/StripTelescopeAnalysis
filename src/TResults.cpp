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
}

TResults::TResults(TSettings *settings) {
  path = gSystem->pwd();
  cout<<"New TResults with settings "<<settings<<"\t run: "<<settings->getRunNumber()<<endl;
  openResults(settings);
}

TResults::~TResults() {
  // TODO Auto-generated destructor stub
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
    stringstream name;
    name << "results_"<<runnumber;
    oldResults = (TResults*)file->FindObject(name.str().c_str());
    cout<<"old Results: "<<oldResults<<endl;
    cout<<oldResults->IsZombie()<<endl;
//    cout<<oldResults->IsA()->ClassName()<<endl;
    oldResults->getLastUpdateDate().Print();
    cout<<"LAST UPDATE ON "<<oldResults->getLastUpdateDate().AsString()<<endl;

  }

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
}
