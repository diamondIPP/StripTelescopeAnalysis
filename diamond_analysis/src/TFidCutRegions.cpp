/*
 * TFidCutRegions.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: bachmair
 */

#include "../include/TFidCutRegions.hh"

using namespace std;
ClassImp(TFidCutRegions);
TFidCutRegions::TFidCutRegions() {
	cout<< "create empty Fid Cut Regions pattern... "<<endl;
	initVariables();
}

void TFidCutRegions::Reset(){
	initVariables();
}
void TFidCutRegions::initVariables(){
	index = 0;
	xInt.clear();
	yInt.clear();
	nFidCuts = 0;
	for(UInt_t i =0;i<fidCuts.size();i++){
		if (fidCuts.at(i))
			delete fidCuts.at(i);
	}
	fidCuts.clear();
	name="";
	verbosity=0;
}

TFidCutRegions::TFidCutRegions(std::vector<std::pair <Float_t, Float_t> > xInt, std::vector<std::pair <Float_t, Float_t> > yInt,UInt_t nDia) {
	initVariables();
	this->nDiamonds=nDia;
	nFidCuts=xInt.size()*yInt.size();
	this->xInt=xInt;
	this->yInt=yInt;
	cout<< "create Fid Cut Regions for a RUN with "<<nDiamonds<<" Diamonds. There are "<<nFidCuts<<" Fiducial Cut Areas"<<endl;
	createFidCuts();
}

TFidCutRegions::TFidCutRegions(Float_t xLow, Float_t xHigh, Float_t yLow, Float_t yHigh, UInt_t nDia){
	initVariables();
	this->nDiamonds=nDia;
	this->addFiducialCut(xLow,xHigh,yLow,yHigh);
}

TFidCutRegions::TFidCutRegions(TH2F* histo,int nDiamonds,Float_t fidCutPercentage)
{
	cout<<"Create FiducialCutRegions for histogram "<<histo->GetTitle()<<" with "<<nDiamonds<<" Diamonds."<<endl;
	initVariables();
	this->nDiamonds=nDiamonds;
	if (histo==0) {
		cerr<<"Histo is not defined"<<endl;
		nFidCuts=0;
		return;
	}
	hEventScatterPlot=histo;
	TH1D* hProjX = hEventScatterPlot->ProjectionX("hFiducialCutOneDiamondHitProjX");
	TH1D* hProjY = hEventScatterPlot->ProjectionY("hFiducialCutOneDiamondHitProjY");
	if(hProjX==0||hProjY==0)
		return;
	hProjX->GetXaxis()->SetTitle("Mean Silicon Value in X[strips]");
	hProjX->GetYaxis()->SetTitle("Number Of Entries #");
	hProjY->GetXaxis()->SetTitle("Mean Silicon Value in Y[strips]");
	hProjY->GetYaxis()->SetTitle("Number Of Entries #");

	hProjX->SetName("hProjX");
	hProjY->SetName("hProjY");
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
	cout<<TCluster::Intent(intend)<<"Printing Fiducial cuts: "<<endl;
	for(UInt_t i=0;i<fidCuts.size();i++){
		fidCuts.at(i)->Print(intend+1);
	}
	cout<<TCluster::Intent(intend)<<"Range to plot all Fiducial Cuts:"<<endl;
	cout<<TCluster::Intent(intend)<<"\tX: "<<getMinFiducialX()<<"-"<<getMaxFiducialX()<<endl;
	cout<<TCluster::Intent(intend)<<"\tY: "<<getMinFiducialY()<<"-"<<getMaxFiducialY()<<endl;

	if(verbosity >3 &&verbosity%2==1){
		cout<<"Done\nPress a key...\t"<<flush;
		char t;
		cin>>t;
	}
}


void TFidCutRegions::setRunDescription(std::string runDes)
{
	index=0;
	this->runDescription=runDes;
	cout<<"runDescription :"<<runDescription<<endl;
	if(runDescription.size()==0) {
		index=0;
		return;
	}
	switch(runDescription.at(0)){
	case '0': index=0;break;//all
	case '1': index=1;break;//left
	case '2': index=2;break;//right
	case '3': index=3;break;//XXX
	case '4': index=4;break;//???
	}

	if(runDescription.find("left")!=string::npos){
		index = 1;
		cout<<" runDescription is :"<<runDescription<<" ==> index: "<<index<<endl;
	}
	else if(runDescription.find("right")!=string::npos){
		index = 2;
		cout<<" runDescription is :"<<runDescription<<" ==> index: "<<index<<endl;
	}
	else if(runDescription.find("middle")!=string::npos){
		index = 3;
		cout<<" runDescription is :"<<runDescription<<" ==> index: "<<index<<endl;
	}

	cout<<" runDescription is :"<<runDescription<<" ==> index: "<<index<<endl;
	if(index!=0){
		if(index>fidCuts.size()){
			cerr<<"the set index ( "<<index<<") is not possible. There are only "<<fidCuts.size()<< " possible Fidcuts..."<<endl;
			cerr<<"Please enter a valid index. Use a number between 0 [all],1(left) - "<<fidCuts.size()+1<<endl;
			unsigned int newIndex=0;
			while(!(cin>>newIndex)||newIndex>fidCuts.size())
			{
				cout<<"you didn't entered a vaild input."<<
						"Please try again. Use a number between 0 [all],1(left) - "<<fidCuts.size()<<endl;
				cin.clear();
				while(cin.get() != '\n'){};
			}
			index = newIndex;
			cout<<" runDescription is :"<<runDescription<<" ==> index: "<<index<<endl;
			return;
		}
		return;
	}
	return;
}

Int_t TFidCutRegions::getFiducialCutIndex(Float_t xVal,Float_t yVal){
	for(UInt_t i=0;i<fidCuts.size();i++)
		if(fidCuts.at(i)->isInFiducialCut(xVal,yVal))
			return i+1;
	return -1;
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
		for(UInt_t i=0;i<fidCuts.size()&&inFidCut==false;i++)
			inFidCut = inFidCut || fidCuts.at(i)->isInFiducialCut(xVal,yVal);
		return inFidCut;
	}
	return false;
}

Float_t TFidCutRegions::getMinFiducialX(UInt_t index){
	if (index ==0){
		Float_t min = 1e37;
		for (UInt_t i = 1;i < this->getNFidCuts()+1; i++)
			min = TMath::Min(min,getMinFiducialX(i));
		return min;
	}
	if(index <= getNFidCuts()){
        //cout<<"GetFidCut: "<<getNFidCuts()<<endl;
		if (getFidCut(index))
			return getFidCut(index)->GetXLow();
	}
	return 1e37;


}

Float_t TFidCutRegions::getMaxFiducialX(UInt_t index){
	if (index == 0){
		Float_t max = -1e37;
		for (UInt_t i = 1;i < this->getNFidCuts()+1; i++)
			max = TMath::Max(max,getMaxFiducialX(i));
		return max;
	}
	if(index <= getNFidCuts()){
		if (getFidCut(index))
			return getFidCut(index)->GetXHigh();
	}
	return -1e37;

}

Float_t TFidCutRegions::getMinFiducialY(UInt_t index){
	if (index ==0){
		Float_t min = 1e37;
		for (UInt_t i = 1;i < this->getNFidCuts()+1; i++)
			min = TMath::Min(min,getMinFiducialY(i));
		return min;
	}
	if(index <= getNFidCuts()){
		if (getFidCut(index))
			return getFidCut(index)->GetYLow();
	}
	return 1e37;

}
Float_t TFidCutRegions::getMaxFiducialY(UInt_t index){
	if (index == 0){
		Float_t max = -1e37;
		for (UInt_t i = 1;i < this->getNFidCuts()+1; i++)
			max = TMath::Max(max,getMaxFiducialY(i));
		return max;
	}
	if(index <= getNFidCuts()){
		if (getFidCut(index))
			return getFidCut(index)->GetYHigh();
	}
	return -1e37;
}
TCanvas *TFidCutRegions::getAllFiducialCutsCanvas(TH2F *hScatter,bool optimizeAxisRange)
{
	if(hScatter==0) hScatter= hEventScatterPlot;
	if(hScatter==0)  return 0;
	TCanvas *c1 = new TCanvas("cFiducialCuts","cFiducialCuts",1024,800);
	c1->cd();
	string hName = hScatter->GetName();
	hName.append("_clone");
	TH1D* projX = hScatter->ProjectionX();
	Float_t binXmin=0;
	for(binXmin=0;binXmin<projX->GetNbinsX()&&projX->GetBinContent(binXmin)==0;binXmin++)
		continue;
	Float_t binXmax=projX->GetNbinsX()-1;
	for(binXmax=projX->GetNbinsX();binXmax>=0&&projX->GetBinContent(binXmax)==0;binXmax--)
		continue;
	if(verbosity){
		cout<<"binXmin = "<<binXmin<<"--->"<<projX->GetBinLowEdge(binXmin)<<endl;
		cout<<"binXmax = "<<binXmax<<"--->"<<projX->GetBinLowEdge(binXmax)<<endl;
	}
	binXmin = projX->GetBinLowEdge(binXmin);
	binXmax = projX->GetBinLowEdge(binXmax);
	delete projX;

	TH1D* projY = hScatter->ProjectionY();
	Float_t binYmin=0;
	for(binYmin=0;binYmin<projY->GetNbinsX()&&projY->GetBinContent(binYmin)==0;binYmin++)
		continue;
	Float_t binYmax=projY->GetNbinsX()-1;
	for(binYmax=projY->GetNbinsX()-1;binYmax>=0&&projY->GetBinContent(binYmax)==0;binYmax--)
		continue;
	if(verbosity){
		cout<<"binYmin = "<<binYmin<<"--->"<<projY->GetBinLowEdge(binYmin)<<endl;
		cout<<"binYmax = "<<binYmax<<"--->"<<projY->GetBinLowEdge(binYmax)<<endl;
	}
	binYmin = projY->GetBinLowEdge(binYmin);
	binYmax = projY->GetBinLowEdge(binYmax);
	delete projY;

	Float_t minX = getMinFiducialX();
	Float_t maxX = getMaxFiducialX();
	Float_t minY = getMinFiducialY();
	Float_t maxY = getMaxFiducialY();
	minX = TMath::Min(minX,binXmin);
	minY = TMath::Min(minY,binYmin);
	maxX = TMath::Max(maxX,binXmax);
	maxY = TMath::Max(maxY,binYmax);
	minX = minX - (maxX-minX)*.1;
	maxX = maxX + (maxX-minX)*.1;
	minY = minY - (maxY-minY)*.1;
	maxY = maxY + (maxY-minY)*.1;
	if(verbosity)
		cout<<"Plot range: x: "<<minX<<"-"<<maxX<<",\ty: "<<minY<<"-"<<maxY<<endl;
	if(verbosity>3&&verbosity%2==1){
		cout<<"\tPress a key and enter to confirm.\t"<<flush;
		char t;
		cin >>t;
	}
	TH2F* htemp = (TH2F*)hScatter->Clone(hName.c_str());
	htemp->Draw("colz");
	if(optimizeAxisRange){
		htemp->GetXaxis()->SetRangeUser(minX,maxX);
		htemp->GetYaxis()->SetRangeUser(minY,maxY);
	}
	for(UInt_t i=0;i<fidCuts.size();i++)
		getFiducialAreaCut(i)->Draw("same");
	return c1;
}

TCutG *TFidCutRegions::getFiducialAreaCut(UInt_t nFidCut)
{
	if (nFidCut>=fidCuts.size())
		return 0;
	Float_t xLow = fidCuts.at(nFidCut)->GetXLow();
	Float_t xHigh = fidCuts.at(nFidCut)->GetXHigh();
	Float_t yLow = fidCuts.at(nFidCut)->GetYLow();
	Float_t yHigh = fidCuts.at(nFidCut)->GetYHigh();
	cout<<"getting fitucial Area for "<<nFidCut<<": ";
	cout<<TString::Format("X: %d-%d, Y: %d-%d",xLow,xHigh,yLow,yHigh);
	TString name = TString::Format("fidCut_%d",nFidCut);
	TCutG * pt = new TCutG(name,5);
	pt->SetPoint(0,xLow,yLow);
	pt->SetPoint(1,xLow,yHigh);
	pt->SetPoint(2,xHigh,yHigh);
	pt->SetPoint(3,xHigh,yLow);
	pt->SetPoint(4,xLow,yLow);
	//(xLow,yLow,xHigh,yHigh);
	if(index == nFidCut + 1||index == 0){
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

int TFidCutRegions::getFidCutRegion(Float_t xVal, Float_t yVal)
{
	unsigned int i=0;
	for(i=0;i<fidCuts.size();i++){
		if(fidCuts.at(i)->isInFiducialCut(xVal,yVal)) return i+1;
	}
	return -1;
}

/**
 * adding a fiducial region
 * @param xLow
 * @param xHigh
 * @param yLow
 * @param yHigh
 * @todo check if region overlaps with another fiducial cut
 */
void TFidCutRegions::addFiducialCut(Float_t xLow, Float_t xHigh, Float_t yLow, Float_t yHigh) {
	if(xLow<=xHigh&&yLow<=yHigh){
		xInt.push_back(make_pair(xLow,xHigh));
		yInt.push_back(make_pair(yLow,yHigh));
		nFidCuts++;

		TFiducialCut* fidCut = new TFiducialCut(nFidCuts,xLow,xHigh,yLow,yHigh);
		fidCuts.push_back(fidCut);
		cout<<name<<": Added another FiducialCut no. "<<nFidCuts<<"_"<<this->getNFidCuts()<<"_"<<fidCuts.size()<<": ["<<xLow<<"-"<<xHigh<<" , "<<yLow<<"-"<<yHigh<<"]"<<endl;
	}
	else
		cerr<<"Could not add Fiducial Cut since "<<xLow<<">="<<xHigh<<" or "<<yLow<<" >= "<<yHigh<<endl;
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

TFiducialCut* TFidCutRegions::getFidCut(UInt_t index){
	if (index>getNFidCuts()&&index>0)
		return NULL;
	if(index==0)
		return NULL;
	return fidCuts.at(index-1);
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
		cout<<"hEventScatterPlot is Zero!!!!!";
		if(verbosity%2==1){
			cout<<"Press any Key to Continue."<<endl;
			char t; cin>>t;
		}
		else
			cout<<endl;
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
	vector<TCutG* > boxes;
	for(UInt_t i=0;i<intervals.size();i++){
		TString name = TString::Format("cut_%s_interval_%d",hProj->GetName(),i);
		TCutG *box = new TCutG(name,5);
		//TPaveText(intervals.at(i).first,0,intervals.at(i).second,hProj->GetMaximum());
		box->SetPoint(0,intervals.at(i).first,0);
		box->SetPoint(1,intervals.at(i).first,hProj->GetMaximum());
		box->SetPoint(2,intervals.at(i).second,hProj->GetMaximum());
		box->SetPoint(3,intervals.at(i).second,0);
		box->SetPoint(4,intervals.at(i).first,0);
		box->SetFillColor(kRed+i);
		box->SetFillStyle(3013);
//		box->AddText("");
//		box->AddText("");
//		box->AddText(Form("Area\n\n%i",i));
		boxes.push_back(box);
		box->Draw("same");
	}
	return c1;

}


Float_t TFidCutRegions::getXLow(UInt_t i) {
	if(i <= getNFidCuts()&&i>0)
		return xInt.at((i-1)%xInt.size()).first;
	return 0;
	return N_INVALID;
}

Float_t TFidCutRegions::getXHigh(UInt_t i) {
	if(i <= getNFidCuts()&&i>0)
		return xInt.at((i-1)%xInt.size()).second;
	return TPlaneProperties::getNChannelsSilicon();
	return N_INVALID;
}

Float_t TFidCutRegions::getYLow(UInt_t i) {
	if(i <= getNFidCuts()&&i>0)
		return xInt.at((i-1)/xInt.size()).first;
	return TPlaneProperties::getNChannelsSilicon();

	return N_INVALID;
}
Float_t TFidCutRegions::getYHigh(UInt_t i) {
	if(i <= getNFidCuts()&&i>0)
		return xInt.at((i-1)/xInt.size()).second;
	return TPlaneProperties::getNChannelsDiamond();
	return N_INVALID;
}
