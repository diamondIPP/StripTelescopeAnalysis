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

	stringstream sectionSummary;
	stringstream section2Highest;
	stringstream sectionAll;
//	sectionSummary<<"<h2>\n"<<"Summary table"<"</h2>\n";
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
//	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) { // DA
	for (UInt_t clusterSize = 0; clusterSize < settings->getMaxTransparentClusterSize(); clusterSize++) { // DA
		vecTable.at(0).push_back(floatToString(clusterSize+1));
		vecTable.at(2).push_back(floatToString(vecMeanPulseHeigths.at(0).at(clusterSize)));
		vecTable.at(3).push_back(floatToString(vecMPPulseHeigths.at(0).at(clusterSize)));
		vecTable.at(5).push_back(floatToString(vecMeanPulseHeigths.at(1).at(clusterSize)));
		vecTable.at(6).push_back(floatToString(vecMPPulseHeigths.at(1).at(clusterSize)));
	}
	sectionSummary << createTable(vecTable);
	sectionSummary	<< putImage(".",(TString)"hDiaTranspAnaPulseHeightMean")
									<< putImage(".",(TString)"hDiaTranspAnaPulseHeightMP")
									<< putImage(".",(TString)"hDiaTranspAnaPulseHeightOf2HighestMean")
									<< putImage(".",(TString)"hDiaTranspAnaPulseHeightOf2HighestMP");
	sectionSummary << "\n\n<br><br>\n\n";
    addSection("Pulse Height - Summary table",sectionSummary.str());
	stringstream plots1, plots2;
//	for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) { // DA
	for (UInt_t clusterSize = 1; clusterSize < settings->getMaxTransparentClusterSize()+1; clusterSize++) { // DA
		stringstream histoname1, histoname2;
		histoname1 << "cDiaTranspAnaPulseHeightOf"<<clusterSize<<"Strips";
		histoname2 << "cDiaTranspAnaPulseHeightOf2HighestIn"<<clusterSize<<"Strips";
		plots1 << putImage(".",histoname1.str()) << " \n";
		plots2 << putImage(".",histoname2.str()) << " \n";
	}
//	sectionAll << "<h2>Pulse Height of N strips</h2><br>";
	sectionAll << plots1.str();
	sectionAll << "\n\n<br><br>\n\n";
//	section2Highest << "<h2>Pulse Height of 2 hightest channels in N strips</h2><br>";
	section2Highest << plots2.str();
	section2Highest << "\n\n<br><br>\n\n";

    addSection("Pulse Height of N strips",sectionAll.str());
    addSection("Pulse Height of 2 hightest channels in N strips",section2Highest.str());

	TString name;
	stringstream pXY;
	name = TString::Format("hLandau2HighestFidCutX_2outOf%02d",settings->getMaxTransparentClusterSize());
//	pXY << putImage((string)".",(string)"hLandau2HighestFidCutX_2outOf10")<<"\n"; // DA
	pXY << putImage((string)".",(string)name) <<"\n"; // DA
	name = TString::Format("hLandau2HighestFidCutY_2outOf%02d",settings->getMaxTransparentClusterSize());
//	pXY << putImage((string)".",(string)"hLandau2HighestFidCutY_2outOf10")<<"\n"; // DA
	pXY << putImage((string)".",(string)name)<<"\n"; // DA
	pXY <<"\n\n<br><br>\n\n";
	name = TString::Format("hLandau2HighestPredHitX_2outOf%02d",settings->getMaxTransparentClusterSize());
//	pXY << putImage((string)".",(string)"hLandau2HighestPredHitX_2outOf10")<<"\n"; // DA
	pXY << putImage((string)".",(string)name)<<"\n"; // DA
	name = TString::Format("hLandau2HighestPredHitY_2outOf%02d",settings->getMaxTransparentClusterSize());
//	pXY << putImage((string)".",(string)"hLandau2HighestPredHitY_2outOf10")<<"\n"; // DA
	pXY << putImage((string)".",(string)name)<<"\n"; // DA
	pXY <<"\n\n<br><br>\n\n";
	name = TString::Format("hLandau2HighestHitProfile_2OutOf%02d",settings->getMaxTransparentClusterSize());
//	pXY << putImage((string)".",(string)"hLandau2HighestHitProfile_2OutOf10")<<"\n"; // DA
	pXY << putImage((string)".",(string)name)<<"\n"; // DA
    addSection("XY Pulse Height Distributions",pXY.str());
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
	vecTable.resize(8);
	vecTable.at(0).push_back("number of used channels");
	vecTable.at(1).push_back("mean using charge weighted position  [&#956m]");
	vecTable.at(2).push_back("resolution [&#956m]");
	vecTable.at(3).push_back("mean using 2 highest channels  [&#956m]");
	vecTable.at(4).push_back("resolution [&#956m]");
	vecTable.at(5).push_back("");
	vecTable.at(6).push_back("mean using 2 highest channels eta corrected  [&#956m]");
	vecTable.at(7).push_back("resolution [&#956m]");
	if(resolutions.size()>3){
		vecTable.resize(10);
		vecTable.at(8).push_back("mean using 2nd Gauss  [&#956m]");
		vecTable.at(9).push_back("resolution [&#956m]");
	}
//	for (UInt_t clusterSize = 0; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) { // DA
	for (UInt_t clusterSize = 0; clusterSize < settings->getMaxTransparentClusterSize(); clusterSize++) { // DA
		vecTable.at(0).push_back(floatToString(clusterSize+1));
		vecTable.at(1).push_back(floatToString(resolutions.at(0).at(clusterSize).first,2));
		vecTable.at(2).push_back(floatToString(resolutions.at(0).at(clusterSize).second,2));
		vecTable.at(3).push_back(floatToString(resolutions.at(1).at(clusterSize).first,2));
		vecTable.at(4).push_back(floatToString(resolutions.at(1).at(clusterSize).second,2));
		vecTable.at(5).push_back("");
		vecTable.at(6).push_back(floatToString(resolutions.at(2).at(clusterSize).first,2));
		vecTable.at(7).push_back(floatToString(resolutions.at(2).at(clusterSize).second,2));
		if(resolutions.size()>3){
			vecTable.at(8).push_back(floatToString(resolutions.at(3).at(clusterSize).first,2));
			vecTable.at(9).push_back(floatToString(resolutions.at(3).at(clusterSize).second,2));
		}
	}
	sectionContent << createTable(vecTable);
	sectionContent << "\n\n<br><br>\n\n";
	stringstream plots1, plots2, plots3;
//	for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector)+1; clusterSize++) { // DA
	for (UInt_t clusterSize = 1; clusterSize < settings->getMaxTransparentClusterSize()+1; clusterSize++) { // DA
		stringstream histoname1, histoname2, histoname3;
		TString name1,name2,name3;
		name1 = TString::Format("hDiaTranspAnaResidualChargeWeightedIn%02dStripsMinusPred",clusterSize);
		name2 = TString::Format("hDiaTranspAnaResidualHighest2CentroidIn%02dStripsMinusPred",clusterSize);
		name3 = TString::Format("hDiaTranspAnaResidualEtaCorrectedIn%02dStripsMinusPred",clusterSize);
		plots1 << putImage(".",(string)name1) << " \n";
		plots2 << putImage(".",(string)name2) << " \n";
		plots3 << putImage(".",(string)name3) << " \n";
	}
	sectionContent << "<h2>Charge weighted position of N strips</h2><br>" << plots1.str();
	sectionContent << "\n\n<br><br>\n\n";
	sectionContent << "<h2>Position of 2 hightest channels in N strips</h2><br>" << plots2.str();
	sectionContent << "\n\n<br><br>\n\n";
	sectionContent << "<h2>Position of eta corrected 2 hightest channels in N strips</h2><br>" << plots3.str();
	addSection("Resolution Plots",sectionContent.str());
}

void THTMLTransparentAnalysis::createEtaPlots() {
	// TODO: change this:
	subjectDetector = 8;

	stringstream sectionContent;
	//	sectionContent<<"<h2>\n"<<
	//	"Eta distri"
	//	<<"</h2>\n";

//	for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) { // DA
	for (UInt_t clusterSize = 1; clusterSize < settings->getMaxTransparentClusterSize(); clusterSize++) { // DA
		stringstream histoname;
		histoname << "hDiaTranspAnaEta2HighestIn"<<clusterSize+1<<"Strips";
		sectionContent << putImage(".",histoname.str());
	}


	addSection("Eta Distributions",sectionContent.str());
}

void THTMLTransparentAnalysis::createEtaIntegrals() {
	// TODO: change this:
	subjectDetector = 8;

	stringstream sectionContent;
	//	sectionContent<<"<h2>\n"<<
	//	"Eta distri"
	//	<<"</h2>\n";

//	for (UInt_t clusterSize = 1; clusterSize < TPlaneProperties::getMaxTransparentClusterSize(subjectDetector); clusterSize++) { // DA
	for (UInt_t clusterSize = 1; clusterSize < settings->getMaxTransparentClusterSize(); clusterSize++) { // DA
		stringstream histoname;
		histoname << "hDiaTranspAnaEtaIntegral2HighestIn"<<clusterSize+1<<"Strips";
		sectionContent << putImage(".",histoname.str());
	}


	addSection("Eta Integrals",sectionContent.str());
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
