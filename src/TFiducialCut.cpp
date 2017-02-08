/*
 * FidCutRegion.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */
#include "TFiducialCut.hh"

ClassImp(TFiducialCut);
using namespace std;
TFiducialCut::TFiducialCut(int i,Float_t xLow,Float_t xHigh,Float_t yLow,Float_t yHigh){
	SetAllValuesZero();
	index=i;
	SetXLow(xLow);
	SetXHigh(xHigh);
	SetYLow(yLow);
	SetYHigh(yHigh);
	name = TString::Format("fidcut_%d",i);
}
TFiducialCut::TFiducialCut(int i) {
	SetAllValuesZero();
	index = i;
	name = TString::Format("fidcut_%d",i);
}

void TFiducialCut::SetAllValuesZero() {
	active = 0;
	x_low = 0;
	x_high = 0;
	y_low = 0;
	y_high = 0;
	name = "fidcut";
}

void TFiducialCut::Print (UInt_t i) {
	std::cout <<TCluster::Intent(i)<< "FidCutRegion #" << index << ":"<<
			" X: " <<TString::Format("%7.1f - %-7.1f", x_low,x_high ) <<
			" Y: " <<TString::Format("%7.1f - %-7.1f", y_low,y_high ) << "\n"<<std::flush;
}

TCutG* TFiducialCut::GetFiducialAreaCut(bool bEmphasis) {
	cout<<"getting fitucial Area for "<<index<<"'"<<name<<"': ";
	cout<<TString::Format("X: %.1f-%.1f, Y: %.1f-%.1f",GetXLow(),GetXHigh(),GetYLow(),GetYHigh());
	TString name = TString::Format("fidCut_%d",index);
	TCutG * pt = new TCutG(name,5);
	pt->SetPoint(0,GetXLow(),GetYLow());
	pt->SetPoint(1,GetXLow(),GetYHigh());
	pt->SetPoint(2,GetXHigh(),GetYHigh());
	pt->SetPoint(3,GetXHigh(),GetYLow());
	pt->SetPoint(4,GetXLow(),GetYLow());
	//(xLow,yLow,xHigh,yHigh);
	if(active||bEmphasis){
		pt->SetFillColor(kRed);
		pt->SetLineColor(kRed);
		pt->SetLineWidth(3);
		pt->SetFillStyle(3013);
	}
	else{
		pt->SetFillColor(kOrange);
		pt->SetFillStyle(3013);
		pt->SetLineColor(kOrange);
		pt->SetLineWidth(3);
	}
	return pt;
}

void TFiducialCut::DrawFiducialCutToCanvas(TCanvas* c1,bool bEmphasis=false) {
	if (!c1)
		return;
	c1->cd();
	TCutG* cut = GetFiducialAreaCut(bEmphasis);
	if (cut)
		cut->Draw("same");

}

void TFiducialCut::SetAllValues(Float_t xLow, Float_t xHigh, Float_t yLow,
        Float_t yHigh) {
    this->SetXLow(xLow);
    this->SetXHigh(xHigh);
    this->SetYLow(yLow);
    this->SetYHigh(yHigh);
}
