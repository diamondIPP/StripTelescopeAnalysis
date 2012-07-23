/*
 * THTMLTransparentAnalysis.cpp
 *
 *  Created on: Jul 5, 2012
 *      Author: bachmair
 */

#include "THTMLTransparentAnalysis.hh"

THTMLTransparentAnalysis::THTMLTransparentAnalysis(TSettings* settings):THTMLGenerator(settings) {
	// TODO Auto-generated constructor stub
	this->setFileName("transparentAnalysis.html");
	this->setMainPath("../");
	this->setSubdirPath("transparentAnalysis/");
	this->setTitle("Transparent Analysis");
	
	
}

THTMLTransparentAnalysis::~THTMLTransparentAnalysis() {
	// TODO Auto-generated destructor stub
}

void THTMLTransparentAnalysis::createContent() {
	
}

void THTMLTransparentAnalysis::createPulseHeightPlots(vector<vector <Float_t> > vecMeanPulseHeigths, vector<vector <Float_t> > vecMPPulseHeigths) {
	// TODO: change this:
	subjectDetector = 8;
	
	stringstream sectionContent;
	sectionContent<<"<h2>\n"<<
	"Summary table"
	<<"</h2>\n";
	std::vector< std::vector< std::string> > vecTable;
//	if(vecMeanPulseHeigths.size()<TPlaneProperties::getNDetectors()) vecMeanPulseHeigths.resize(TPlaneProperties::getNDetectors());
	vecTable.resize(7);
	vecTable.at(0).push_back("number of used channels");
	vecTable.at(1).push_back("PulseHeigth");
	vecTable.at(2).push_back("mean");
	vecTable.at(3).push_back("most probable");
	vecTable.at(4).push_back("PulseHeigth 2 highest channels");
	vecTable.at(5).push_back("mean");
	vecTable.at(6).push_back("most probable");
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		vecTable.at(0).push_back(floatToString(clusterSize+1));
		vecTable.at(2).push_back(floatToString(vecMeanPulseHeigths.at(0).at(clusterSize)));
		vecTable.at(3).push_back(floatToString(vecMPPulseHeigths.at(0).at(clusterSize)));
		vecTable.at(5).push_back(floatToString(vecMeanPulseHeigths.at(1).at(clusterSize)));
		vecTable.at(6).push_back(floatToString(vecMPPulseHeigths.at(1).at(clusterSize)));
	}
	sectionContent << createTable(vecTable);
	sectionContent	<< putImage(".","hDiaTranspAnaPulseHeightMean")
					<< putImage(".","hDiaTranspAnaPulseHeightMP")
					<< putImage(".","hDiaTranspAnaPulseHeightOf2HighestMean")
					<< putImage(".","hDiaTranspAnaPulseHeightOf2HighestMP");
	sectionContent << "\n\n<br><br>\n\n";
	stringstream plots1, plots2;
	for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) {
		stringstream histoname1, histoname2;
		histoname1 << "hDiaTranspAnaPulseHeightOf"<<clusterSize<<"Strips";
		histoname2 << "hDiaTranspAnaPulseHeightOf2HighestIn"<<clusterSize<<"Strips";
		plots1 << putImage(".",histoname1.str()) << " \n";
		plots2 << putImage(".",histoname2.str()) << " \n";
	}
	sectionContent << "<h2>Pulse Height of N strips</h2><br>" << plots1.str();
	sectionContent << "\n\n<br><br>\n\n";
	sectionContent << "<h2>Pulse Height of 2 hightest channels in N strips</h2><br>" << plots2.str();
	addSection("Pulse Height Distributions",sectionContent.str());
}

void THTMLTransparentAnalysis::createResolutionPlots(vector<vector <pair <Float_t,Float_t> > > resolutions) {
	// TODO: change this:
	subjectDetector = 8;
	
	stringstream sectionContent;
	sectionContent<<"<h2>\n"<<
	"Summary table"
	<<"</h2>\n";
	std::vector< std::vector< std::string> > vecTable;
	//	if(vecMeanPulseHeigths.size()<TPlaneProperties::getNDetectors()) vecMeanPulseHeigths.resize(TPlaneProperties::getNDetectors());
	vecTable.resize(5);
	vecTable.at(0).push_back("number of used channels");
	vecTable.at(1).push_back("mean & width using charge weighted position");
	vecTable.at(2).push_back("resolution [&#956m]");
	vecTable.at(3).push_back("mean & width using 2 highest channels");
	vecTable.at(4).push_back("resolution [&#956m]");
	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) {
		vecTable.at(0).push_back(floatToString(clusterSize+1));
		vecTable.at(1).push_back(floatToString(resolutions.at(0).at(clusterSize).first)+"|"+floatToString(resolutions.at(0).at(clusterSize).second));
		vecTable.at(2).push_back(floatToString(50.*resolutions.at(0).at(clusterSize).second));
		vecTable.at(3).push_back(floatToString(resolutions.at(1).at(clusterSize).first)+"|"+floatToString(resolutions.at(1).at(clusterSize).second));
		vecTable.at(4).push_back(floatToString(50.*resolutions.at(1).at(clusterSize).second));
	}
	sectionContent << createTable(vecTable);
	sectionContent << "\n\n<br><br>\n\n";
	stringstream plots1, plots2;
	for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) {
		stringstream histoname1, histoname2;
		histoname1 << "hDiaTranspAnaResidualChargeWeightedIn"<<clusterSize<<"StripsMinusPred";
		histoname2 << "hDiaTranspAnaResidualHighest2CentroidIn"<<clusterSize<<"StripsMinusPred";
		plots1 << putImage(".",histoname1.str()) << " \n";
		plots2 << putImage(".",histoname2.str()) << " \n";
	}
	sectionContent << "<h2>Charge weighted position of N strips</h2><br>" << plots1.str();
	sectionContent << "\n\n<br><br>\n\n";
	sectionContent << "<h2>Position of 2 hightest channels in N strips</h2><br>" << plots2.str();
	addSection("Resolution Plots",sectionContent.str());
}

void THTMLTransparentAnalysis::createEtaPlots() {
	
}


//std::string THTMLGenerator::putImagesOfAllDetectors(std::string path,std::string name, std::string type,int percentage){
//	
//	stringstream output;
//	output<<"\n\t";
//	for(UInt_t det = 0; det< TPlaneProperties::getNSiliconDetectors();det+=2){
//		stringstream name2;
//		name2<<name<<TPlaneProperties::getStringForDetector(det);
//		output<<putImage(path,name2.str());
//	}
//	output<<"\n<br\n\t";
//	for(UInt_t det = 1; det< TPlaneProperties::getNSiliconDetectors();det+=2){
//		stringstream name2;
//		name2<<name<<TPlaneProperties::getStringForDetector(det);
//		output<<putImage(path,name2.str());
//	}
//	output<<"\n<br>\n\t";
//	stringstream name2;
//	name2<<name<<TPlaneProperties::getStringForDetector(TPlaneProperties::getDetDiamond());
//	output<<putImage(path,name2.str());
//	return (output.str());
//}