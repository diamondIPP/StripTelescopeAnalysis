/*
 * TFidCutRegions.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: bachmair
 */

#include "../include/TFidCutRegions.hh"

using namespace std;
TFidCutRegions::TFidCutRegions(std::vector<std::pair <Float_t, Float_t> > xInt,std::vector<std::pair <Float_t, Float_t> >yInt,UInt_t nDia) {
  this->nDiamonds=nDia;
  index=0;
  nFidCuts=xInt.size()*yInt.size();
  this->xInt=xInt;
  this->yInt=yInt;

  cout<< "create Fid Cut Regions for a RUN with "<<nDiamonds<<" Diamonds. There are "<<nFidCuts<<" Fiducial Cut Areas"<<endl;
  createFidCuts();
}

TFidCutRegions::TFidCutRegions(TH1D* histo, int nDiamonds)
{
  index=0;
  this->nDiamonds=nDiamonds;
//todo create a way that this function extract the fidCuts out of the histo
}


TFidCutRegions::~TFidCutRegions() {
  for(UInt_t i=0;i<fidCuts.size();i++)
    delete fidCuts.at(i);
  fidCuts.clear();
}


void TFidCutRegions::Print(int intend)
{
  cout<<"Printing Fiducial cuts: "<<endl;
  for(UInt_t i=0;i<fidCuts.size();i++){
    fidCuts.at(i)->Print();
  }
}


void TFidCutRegions::setRunDescription(std::string runDes)
{
  this->runDescription=runDes;
  if(runDescription.size()==0) {
    index=0;
    return;
  }
  switch(runDescription.at(0)){
    case '0': index=0;return;
    case '1': index=1;return;
    case '2': index=2;return;
    case '3': index=3;return;
    case '4': index=4;return;
  }
  if(runDescription.find("left")!=string::npos){
    index = 1;
    return;
  }
  if(runDescription.find("right")!=string::npos){
    index = 2;
    return;
  }
  return;
}

bool TFidCutRegions::isInFiducialCut(Float_t xVal, Float_t yVal)
{
//  cout<<"xVal = "<<xVal<<"\tyVal = "<<yVal<<"\tIndex = "<<index<<endl;
  if(index>0){
    if( index<this->fidCuts.size()+1){
//      TFiducialCut fidCut = fidCuts.at(index-1);
      return fidCuts.at(index-1)->isInFiducialCut(xVal,yVal);
    }
  }
  if(index==0){
    bool inFidCut=false;
    for(UInt_t i=0;i<fidCuts.size();i++)
      inFidCut = inFidCut || fidCuts.at(i)->isInFiducialCut(xVal,yVal);
    return inFidCut;
  }
  return false;
}

void TFidCutRegions::createFidCuts(){
  if(nDiamonds!=nFidCuts){
    cout<<"Fid Cut does not match with nDIamonds"<<endl;
  }
  cout<<"\ncreate FidCuts"<<endl;
  for(UInt_t iY=0;iY<yInt.size();iY++)
    for(UInt_t iX=0;iX<xInt.size();iX++){
      cout<<"iX="<<iX<<"\t"<<"iY"<<iY<<" - "<<flush;
      int i = iY*xInt.size()+iX;
      int xLow=xInt.at(iX).first;
      int xHigh=xInt.at(iX).second;
      int yLow = yInt.at(iY).first;
      int yHigh = yInt.at(iY).second;
      TFiducialCut* fidCut = new TFiducialCut(i,xLow,xHigh,yLow,yHigh);
      fidCut->Print();
      this->fidCuts.push_back(fidCut);
    }
}

TFiducialCut* TFidCutRegions::getFidCut(std::string describtion){
  if(fidCuts.size()==0){
    cout<<"fidCuts not yet defined..."<<endl;
    return 0;
  }
  if(describtion.at(0)=='0'){
    cout<<"return empty "<<endl;
    return fidCuts.at(0);
  }
  else if(describtion.find("left")!=string::npos||describtion.at(0)=='1'){
      cout<<"FidCut return at 0"<<endl;
      return fidCuts.at(0);
  }
  else if(describtion.find("right")!=string::npos||describtion.at(0)=='0'){
    if(fidCuts.size()>1)
    return fidCuts.at(1);
    else return fidCuts.at(0);
  }
  return 0;


}
