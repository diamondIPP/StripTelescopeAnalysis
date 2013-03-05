/*
 * THTMLLandaus.cpp
 *
 *  Created on: May 20, 2012
 *      Author: bachmair
 */

#include "../include/THTMLLandaus.hh"

THTMLLandaus::THTMLLandaus(TSettings* settings):THTMLGenerator(settings) {
	this->setFileName("landaus.html");
	this->setMainPath("../");
	this->setSubdirPath("selectionAnalysis/");
	this->setTitle("Landau Distributions");

}

THTMLLandaus::~THTMLLandaus() {
}


void THTMLLandaus::addLandauDiamondTable(vector<Float_t> vecHistoMeans, vector<Float_t> vecHistoMaxs, vector<Float_t> vecHistoGaus, vector<Float_t> vecHistoLandau)
{
	stringstream sectionContent;
	vector<vector<string> > tableContent;
	tableContent.resize(vecHistoMeans.size()+1);
	tableContent.at(0).push_back("ClusterSize");
	tableContent.at(0).push_back("Mean");
	tableContent.at(0).push_back("MaxPos");
	tableContent.at(0).push_back("GausPos");
	tableContent.at(0).push_back("LandauMP");
	tableContent.at(0).push_back("Landau/Gaus");
	tableContent.at(0).push_back("Landau/Max");
	for(UInt_t i=0;i<vecHistoMeans.size()&&i<vecHistoMaxs.size()&&i<vecHistoGaus.size()&&i<vecHistoLandau.size();i++){
		Float_t fraction = vecHistoLandau.at(i)/vecHistoGaus.at(i);
		tableContent.at(i+1).push_back((i==0?"AllClusters":this->floatToString(i+1)));
		tableContent.at(i+1).push_back(this->floatToString(vecHistoMeans.at(i)));
		tableContent.at(i+1).push_back(this->floatToString(vecHistoMaxs.at(i)));
		tableContent.at(i+1).push_back(this->floatToString(vecHistoGaus.at(i)));
		tableContent.at(i+1).push_back(this->floatToString(vecHistoLandau.at(i)));
		tableContent.at(i+1).push_back(this->floatToString(fraction));
		fraction = vecHistoLandau.at(i)/vecHistoMaxs.at(i);
		tableContent.at(i+1).push_back(this->floatToString(fraction));
	}
	sectionContent<<"<p>\n";
	sectionContent<<createTable(tableContent)<<endl;
	sectionContent<<"<p>";
	this->addSection("Landau Fit CrossCheck Table",sectionContent.str());
}

void THTMLLandaus::addLandauDiamond(Float_t width, Float_t MP, Float_t area, Float_t GSigma)
{
	stringstream sectionContent;
	sectionContent<<"<p>\n";
	sectionContent<<"Two dimensional distribution of ClusterCharge vs. ClusterSize for all Events which are used\n";
	sectionContent<<"for Alignment or Analysis, i.e. one and only one cluster in ALL detectors<br>\n";
	sectionContent<<putImage(this->path,(TString)"hLandauDiamond_OneCluster","png",50)<<"<br>\n";
	sectionContent<<"<br>\nThe next plot shows the Projection of the two dimensional plot to the x axis.<br\n";
	sectionContent<<"This is the Plot for all ClusterSizes. ";
	sectionContent<<putImage(this->path,(TString)"hPulseHeightDiamondAll","png",49)<<" \n";
	sectionContent<<putImage(this->path,(TString)"hPulseHeigthDiamond_1_2_ClusterSize","png",49)<<"<br>\n";
	for(UInt_t i=1;i<8;i++){
		stringstream name;
		name <<"hPulseHeigthDiamond_"<<i<<"_ClusterSize";
		sectionContent<<putImage(this->path,name.str(),"png",24)<<"\n"<<(i%4==0?"<br>":"");
	}
	sectionContent<<putImage(this->path,(TString)"cMVP_Landau_vs_ClusterSize","png",50)<<"<br>\n";
	//	sectionContent<<putImage(this->path,"hLandauDiamond_OneCluster","png",50)<<"<br>\n";
	sectionContent<<"</p>";

	this->addSection("Landau Distributions Diamond",sectionContent.str());

}



