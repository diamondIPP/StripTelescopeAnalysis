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

TFidCutRegions::TFidCutRegions(TH2F* histo,int nDiamonds,Float_t fidCutPercentage)
{
  cout<<"Create FiducialCutRegions for histogram "<<histo->GetTitle()<<" with "<<nDiamonds<<" Diamonds."<<endl;
  index=0;
  this->nDiamonds=nDiamonds;
  if (histo==0) {
    cerr<<"Histo is not defined"<<endl;
    nFidCuts=0;
    return;
  }
  hEventScatterPlot=histo;
  TH1D* hProjX = hEventScatterPlot->ProjectionX("hFiducialCutOneDiamondHitProjX");
  hProjX->GetXaxis()->SetTitle("Mean Silicon Value in X[strips]");
  hProjX->GetYaxis()->SetTitle("Number Of Entries #");
  TH1D* hProjY = hEventScatterPlot->ProjectionY("hFiducialCutOneDiamondHitProjY");
  hProjY->GetXaxis()->SetTitle("Mean Silicon Value in Y[strips]");
  hProjY->GetYaxis()->SetTitle("Number Of Entries #");

  hProjX->SetName("hProjX");
  hProjY->SetName("hProjY");
  if(hProjX==0||hProjY==0)
    return;
  if(hProjX->IsZombie()||hProjY->IsZombie())
    return;

  vector<pair <Float_t, Float_t> >xIntervals = findFiducialCutIntervall( hProjX,fidCutPercentage);
  vector<pair <Float_t, Float_t> >yIntervals = findFiducialCutIntervall( hProjY,fidCutPercentage);
  xInt = xIntervals;
  yInt = yIntervals;
  createFidCuts();
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
	index=0;
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
  if(index>0){// one special diamond is choosen e.g. left/right
    if( index<this->fidCuts.size()+1){
//      TFiducialCut fidCut = fidCuts.at(index-1);
//      cout<<xVal<<"/"<<yVal<<"\t"<<index<<" "<<fidCuts.at(index-1)->isInFiducialCut(xVal,yVal)<<endl;
      return fidCuts.at(index-1)->isInFiducialCut(xVal,yVal);
    }
  }
  if(index==0){ //all diamonds will be analyesd
    bool inFidCut=false;
    for(UInt_t i=0;i<fidCuts.size();i++)
      inFidCut = inFidCut || fidCuts.at(i)->isInFiducialCut(xVal,yVal);
    return inFidCut;
  }
  return false;
}

TCanvas *TFidCutRegions::getAllFiducialCutsCanvas(TH2F *hScatter)
{
  if(hScatter==0) hScatter= hEventScatterPlot;
  if(hScatter==0)  return 0;
  TCanvas *c1 = new TCanvas("cFiducialCuts","cFiducialCuts",1024,800);
  c1->cd();
  string hName = hScatter->GetName();
  hName.append("_clone");
  TH2F* htemp = (TH2F*)hScatter->Clone(hName.c_str());
  htemp->Draw("colz");
  for(UInt_t i=0;i<fidCuts.size();i++)
    getFiducialAreaPaveText(i)->Draw("same");
  return c1;
}

TPaveText *TFidCutRegions::getFiducialAreaPaveText(UInt_t nFidCut)
{
  if (nFidCut>=fidCuts.size())
    return 0;
  Float_t xLow = fidCuts.at(nFidCut)->GetXLow();
  Float_t xHigh = fidCuts.at(nFidCut)->GetXHigh();
  Float_t yLow = fidCuts.at(nFidCut)->GetYLow();
  Float_t yHigh = fidCuts.at(nFidCut)->GetYHigh();
  TPaveText * pt = new TPaveText(xLow,yLow,xHigh,yHigh);
  if(index == nFidCut + 1||index == 0){
    pt->SetFillColor(kRed);
    pt->SetFillStyle(3013);
  }
  else{
    pt->SetFillColor(kOrange);
    pt->SetFillStyle(3013);
  }
  return pt;
}

int TFidCutRegions::getFidCutRegion(Float_t xVal, Float_t yVal)
{
  int i=0;
  for(i=0;i<fidCuts.size();i++){
    if(fidCuts.at(i)->isInFiducialCut(xVal,yVal)) return i;
  }
  return -1;
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
  cout<<"DONE with creating FidCuts"<<endl;
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


std::vector< std::pair< Float_t,Float_t> > TFidCutRegions::findFiducialCutIntervall(TH1D* hProj,Float_t fidCutPercentage){
  int minWidth = 03;//channels
  Float_t mean = hProj->Integral()/hProj->GetNbinsX();
  Int_t nLow=0;
  Int_t nHigh=0;
  int nbins = hProj->GetNbinsX();
  Float_t value = mean*fidCutPercentage;
  cout<<value<<endl;
  std::vector< std::pair<Float_t,Float_t> > intervals;
  while(nLow<nbins-1){
    nLow=nHigh+1;
    bool foundLeft=false;
    bool foundRight=false;
    while(nLow<nbins-1&&hProj->GetBinContent(nLow)<value){
      nLow++;
    }
    nHigh=nLow;
    foundLeft = (hProj->GetBinContent(nLow)>value&&nLow<nbins);
    if(foundLeft){
      while(nHigh<nbins&&hProj->GetBinContent(nHigh)>value){
        nHigh++;
      }
      foundRight=hProj->GetBinContent(nHigh)<value&&nHigh<nbins-1;
    }

    if(foundLeft&&foundRight){
      Float_t x1=hProj->GetBinLowEdge(nLow+1)+1;
      Float_t x2=hProj->GetBinLowEdge(nHigh-1)-1;
      if(x1<x2&&nHigh-nLow>minWidth){
        cout<<"Found first intervall: ["<<nLow<<","<<nHigh<<"]\t= ["<<x1<<","<<x2<<"]"<<endl;
        intervals.push_back( std::make_pair(x1,x2));
      }
    }
    nLow=nHigh;
  }
//  DrawFiduciaCuts(hProj,intervals);
//  cout<<"DONE"<<endl;
  return intervals;
}


TCanvas *TFidCutRegions::getFiducialCutCanvas(TPlaneProperties::enumCoordinate cor){
  if(hEventScatterPlot==0) {
    cout<<"hEventScatterPlot is Zero!!!!!"<<endl;
    char t; cin>>t;
    return 0;
  }
  if(hEventScatterPlot->IsZombie())return 0;
  if(cor == TPlaneProperties::X_COR){
   TH1D *hProjX =  hEventScatterPlot->ProjectionX("hFiducialCutOneDiamondHitProjX");
   hProjX->GetXaxis()->SetTitle("Mean Silicon Value in X[strips]");
   hProjX->GetYaxis()->SetTitle("Number Of Entries #");
   TCanvas *c1 = getFiducialCutProjectionCanvas(hProjX,xInt);
   c1->SetName("chProjX");
   return c1;
  }
  else if(cor == TPlaneProperties::Y_COR){
    TH1D *hProjY = hEventScatterPlot->ProjectionY("hFiducialCutOneDiamondHitProjY");
    hProjY->GetXaxis()->SetTitle("Mean Silicon Value in Y[strips]");
    hProjY->GetYaxis()->SetTitle("Number Of Entries #");
    TCanvas *c1 = getFiducialCutProjectionCanvas(hProjY,yInt);
       c1->SetName("chProjY");
       return c1;
  }
  else if(cor == TPlaneProperties::XY_COR){
    return getAllFiducialCutsCanvas();
  }
  return 0;
}

TCanvas *TFidCutRegions::getFiducialCutProjectionCanvas(TH1D* hProj,std::vector< std::pair<Float_t,Float_t> > intervals){
  if (hProj==0)
    return 0;
  std::stringstream canvasName;
  canvasName<<"c"<<hProj->GetName();
  TCanvas *c1 =new TCanvas(canvasName.str().c_str(),canvasName.str().c_str(),800,600);
  c1->cd();
  hProj->Draw();
  vector<TBox* > boxes;
  for(UInt_t i=0;i<intervals.size();i++){
    TPaveText *box = new TPaveText(intervals.at(i).first,0,intervals.at(i).second,hProj->GetMaximum());
    box->SetFillColor(kRed+i);
    box->SetFillStyle(3013);
    box->AddText("");
    box->AddText("");
    box->AddText(Form("Area\n\n%i",i));
    boxes.push_back(box);
    box->Draw("same");
  }
  return c1;

}
