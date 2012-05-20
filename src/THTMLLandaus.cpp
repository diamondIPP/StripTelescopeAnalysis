/*
 * THTMLLandaus.cpp
 *
 *  Created on: May 20, 2012
 *      Author: bachmair
 */

#include "../include/THTMLLandaus.hh"

THTMLLandaus::THTMLLandaus(TSettings* settings):THTMLGenerator(settings) {
	// TODO Auto-generated constructor stub

	this->setFileName("landaus.html");
	this->setSubdirPath("/selectionAnalysis/");
	this->setTitle("Landau Distributions");

}

THTMLLandaus::~THTMLLandaus() {
	// TODO Auto-generated destructor stub
}

void THTMLLandaus::addLandauDiamond(Float_t width, Float_t MP, Float_t area, Float_t GSigma)
{
	stringstream sectionContent;
	sectionContent<<"<p>";
	sectionContent<<putImage(this->path,"hLandauDiamond_OneCluster","png",50)<<"<br>\n";
	sectionContent<<putImage(this->path,"hPulseHeightDiamondAll","png",50)<<"<br>\n";
	for(UInt_t i=1;i<8;i++){
		stringstream name;
		name <<"hPulseHeigthDiamond_"<<i<<"_ClusterSize";
		sectionContent<<putImage(this->path,name.str(),"png",24)<<"\n"<<(i%4==0?"<br>":"");
	}
	sectionContent<<putImage(this->path,"hMPV_Landau_diff_ClusterSizes","png",50)<<"<br>\n";
//	sectionContent<<putImage(this->path,"hLandauDiamond_OneCluster","png",50)<<"<br>\n";
	sectionContent<<"</p>";

	this->addSection("Landau Distributions Diamond",sectionContent.str());

}



