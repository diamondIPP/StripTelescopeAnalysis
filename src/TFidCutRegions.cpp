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
  nFidCuts=xInt.size()*yInt.size();
  this->xInt=xInt;
  this->yInt=yInt;
  cout<< "create Fid Cut Regions for a RUN with "<<nDiamonds<<" Diamonds. There are "<<nFidCuts<<" Fiducial Cut Areas";
  createFidCuts();
}

TFidCutRegions::TFidCutRegions(TH1D* histo, int nDiamonds)
{
  this->nDiamonds=nDiamonds;
//todo create a way that this function extract
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


void TFidCutRegions::createFidCuts(){
  if(nDiamonds==nFidCuts){
    cout<<"\ncreate FidCuts"<<endl;
    for(UInt_t iY=0;iY<yInt.size();iY++)
    for(UInt_t iX=0;iX<xInt.size();iX++){
      cout<<"iX="<<iX<<"\t"<<"iY"<<iY<<" - "<<flush;
      int i = iY*xInt.size()+iX;
      int xLow=xInt.at(iX).first;
      int xHigh=xInt.at(iX).second;
      int yLow = yInt.at(iY).first;
      int yHigh = yInt.at(iY).second;
      TFiducialCut* fidCut = new TFiducialCut(i,xLow,xHigh,yLow,yHigh);//iY*xInt.size()+iX,xInt.at(iX).first,xInt.at(iX).second,yInt.at(iY).first,yInt.at(iY).second);
      cout<<"$ "<<endl;
      fidCut->Print();
      this->fidCuts.push_back(fidCut);
    }
  }
}

TFiducialCut* TFidCutRegions::getFidCut(std::string describtion){
  if(fidCuts.size()==0){
    cout<<"fidCuts not yet defined..."<<endl;
  }
  if(describtion.at(0)=='0'){
    return 0;
  }
  else if(describtion.find("left")!=string::npos||describtion.at(0)=='1'){
      return fidCuts.at(0);
  }
  else if(describtion.find("right")!=string::npos||describtion.at(0)=='0'){
    return fidCuts.at(1);
  }
  return 0;


}
