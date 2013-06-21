/**
 * @file TCellAnalysisClass.cpp
 *
 * @date Jun 20, 2013
 * @author bachmair
 * @description
 */

#include "../include/TCellAnalysisClass.hh"

TCellAnalysisClass::TCellAnalysisClass(TSettings * set) {
	if(set!=0)
		settings=set;
	else
		exit(-1);
	pCluster=0;
	// TODO Auto-generated constructor stub
	Float_t nan = 1./0.;
	cellAnalysisTree = new TTree("cellAnalysisTree","CellAnalysisTree");
	initialiseBranchAddresses();

}

TCellAnalysisClass::~TCellAnalysisClass() {
	// TODO Auto-generated destructor stub
}

void TCellAnalysisClass::addEvent(Float_t xPred, Float_t yPred, Int_t cell,
		UInt_t quarter, Float_t relCellPosX, Float_t relCellPosY,TCluster clus) {
	ph = clus.getCharge();
	clusterSize = clus.getClusterSize();
	this->xPred = xPred;
	this->yPred = yPred;
	this->relCellPosX = relCellPosX;
	this->relCellPosY = relCellPosY;
	this->nQuarter = quarter;
	this->nCell = cell;
	this->nColumn = cell/settings->getNRows3d();
	this->nRow = cell % settings->getNRows3d();
	if(nColumn<0|| nColumn>12){
		cout<<"cannot fill column/row: "<<cell<<" = "<<nColumn<<"/"<<nRow<<endl;
		nColumn = -1;
	}
	if(nRow<0||nRow>12){
		cout<<"cannot fill column/row: "<<cell<<" = "<<nColumn<<"/"<<nRow<<endl;
		nRow = -1;
	}
	cellAnalysisTree->Fill();
}

TH1* TCellAnalysisClass::getHistogram(string histoName,string varexp, string selection, string drawOption) {
//	TH1* histo = new TH1("histo","histo");
	varexp.append(">>");
	varexp.append(histoName);
	cellAnalysisTree->Draw(varexp.c_str(),selection.c_str());
	TH1* histo = (TH1*) gROOT->FindObject(histoName.c_str());
	return histo;
}

void TCellAnalysisClass::initialiseBranchAddresses() {
	if(!cellAnalysisTree)
		return;
//	cellAnalysisTree->Branch("D0X_ADC",&Det_ADC[0],"D0X_ADC[256]/b");
//	cellAnalysisTree->Branch("pulseHeight", "Event", &event, bsize,split);
	cellAnalysisTree->Branch("pulseHeight",&ph,"pulseHeight/F");
	cellAnalysisTree->Branch("xPred",&xPred,"xPred/F");
	cellAnalysisTree->Branch("yPred",&yPred,"yPred/F");
	cellAnalysisTree->Branch("nCell",&nCell,"nCell/I");
	cellAnalysisTree->Branch("nRow",&nRow,"nRow/I");
	cellAnalysisTree->Branch("nColumn",&nColumn,"nColumn/I");
	cellAnalysisTree->Branch("nQuarter",&nQuarter,"nQuarter/I");
	cellAnalysisTree->Branch("relCellPosX",&relCellPosX,"relCellPosX/F");
	cellAnalysisTree->Branch("relCellPosY",&relCellPosY,"relCellPosY/F");
	cellAnalysisTree->Branch("clusterSize",&clusterSize,"clusterSize/I");
	pCluster=0;
	cellAnalysisTree->Branch("diaCluster","TCluster",&pCluster);
}
