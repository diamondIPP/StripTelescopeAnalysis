/*
 * THTMLPedestal.cpp
 *
 *  Created on: May 2, 2012
 *      Author: bachmair
 */

#include "../include/THTMLPedestal.hh"

THTMLPedestal::THTMLPedestal(TSettings *settings):THTMLGenerator(settings) {
	// TODO Auto-generated constructor stub
	setTitle("Pedestals");
}

THTMLPedestal::~THTMLPedestal() {
	// TODO Auto-generated destructor stub
}

void THTMLPedestal::createTableOfCuts()
{
	stringstream sectionContent;
	sectionContent<<"<h2>Seed and Hit Values in Units of Sigma</h2>\n";
	std::vector<std::vector< std::string > > tablecontent;
	std::vector<std::vector< std::string > > tablecontent2;
	tablecontent.resize(3);
	tablecontent.at(0).push_back("Detector");
	tablecontent.at(1).push_back("Seed");
	tablecontent.at(2).push_back("Hit");
	tablecontent2.resize(3);
	tablecontent2.at(0).push_back("Detector");
	tablecontent2.at(1).push_back("Seed");
	tablecontent2.at(2).push_back("Hit");
	for(UInt_t det =0;det <TPlaneProperties::getNDetectors();det+=2){
		tablecontent.at(0).push_back(TADCEventReader::getStringForDetector(det));
		tablecontent.at(1).push_back(floatToString(settings->getClusterSeedFactor(det)));
		tablecontent.at(2).push_back(floatToString(settings->getClusterHitFactor(det)));
	}
	for(UInt_t det =1;det <TPlaneProperties::getNDetectors();det+=2){
		tablecontent2.at(0).push_back(TADCEventReader::getStringForDetector(det));
		tablecontent2.at(1).push_back(floatToString(settings->getClusterSeedFactor(det)));
		tablecontent2.at(2).push_back(floatToString(settings->getClusterHitFactor(det)));
	}
	sectionContent<<"<br><h4> X Coordinates</h4><br>"<<this->createTable(tablecontent)<<"<br><br>";
	sectionContent<<"<br><h4> Y Coordinates</h4><br>"<<this->createTable(tablecontent2)<<"<br><br>";
	stringstream path;
	path<<this->path<<"/";//"./pedestalAnalysis/";
	sectionContent<<"<h3>Seed-Cuts</h3>\n";
	sectionContent<<putImagesOfAllDetectors(path.str(),"hPulseHeight_BiggestHitChannelInSigma");
//	for(UInt_t det = 0; det< TPlaneProperties::getNSiliconDetectors();det+=2){
//		stringstream name;
//		name<<"hPulseHeight_BiggestHitChannelInSigma"<<TADCEventReader::getStringForDetector(det);
//		sectionContent<<putImage(path.str(),name.str());
//	}
//	for(UInt_t det = 1; det< TPlaneProperties::getNSiliconDetectors();det+=2){
//		stringstream name;
//		name<<"hPulseHeight_BiggestHitChannelInSigma"<<TADCEventReader::getStringForDetector(det);
//		sectionContent<<putImage(path.str(),name.str());
//	}
//	stringstream name;
//	name<<"hPulseHeight_BiggestHitChannelInSigma"<<TADCEventReader::getStringForDetector(TPlaneProperties::getDetDiamond());
//	sectionContent<<putImage(path.str(),name.str());
	sectionContent<<"<h3>Hit-Cuts</h3>\n";
	sectionContent<<putImagesOfAllDetectors(path.str(),"hPulseHeight_SecondBiggestHitChannelInSigma_");
//	for(UInt_t det = 0; det< TPlaneProperties::getNSiliconDetectors();det+=2){
//		stringstream name;
//		name<<"hPulseHeight_SecondBiggestHitChannelInSigma_"<<TADCEventReader::getStringForDetector(det);
//		sectionContent<<putImage(path.str(),name.str());
//	}
//	sectionContent<<"<br";
//	for(UInt_t det = 1; det< TPlaneProperties::getNSiliconDetectors();det+=2){
//		stringstream name;
//		name<<"hPulseHeight_SecondBiggestHitChannelInSigma_"<<TADCEventReader::getStringForDetector(det);
//		sectionContent<<putImage(path.str(),name.str());
//	}
//	name.str("");name.clear();name.str("");
//	name<<"hPulseHeight_SecondBiggestHitChannelInSigma_"<<TADCEventReader::getStringForDetector(TPlaneProperties::getDetDiamond());
//	sectionContent<<putImage(path.str(),name.str());
	this->addSection("Cluster Cuts",sectionContent.str());
}

void THTMLPedestal::createPedestalDistribution(){

	stringstream sectionContent;
	sectionContent<<"<h2> Pedestal Distribution</h2>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"Mean pedestal value of every channel calculated for all events in black.";
	sectionContent<<"The mean pedestal sigma of every channel is plotted in red.\n";
	sectionContent<<"</p>\n";
	stringstream path;
	path<<this->path<<"/";
	sectionContent<<putImagesOfAllDetectors(path.str(),"cPedestalOfChannels_");

	this->addSection("mean Pedestal Values",sectionContent.str());

}

void THTMLPedestal::createPageContent()
{
	createTableOfCuts();
	createPedestalDistribution();
	createBiggestHitMaps();
	createNoiseDistribution();
	createHitOrderSection();
	createSaturatedChannels();
}

void THTMLPedestal::createBiggestHitMaps()
{

	stringstream sectionContent;
	sectionContent<<"<h2> Biggest Hit Maps</h2>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"Channel position of Biggest Hit in detector of each event.\n";
	sectionContent<<"</p>\n";
	stringstream path;
	path<<this->path<<"/";
	sectionContent<<putImagesOfAllDetectors(path.str(),"hBiggestHitMap");
	this->addSection("Biggest Hit Maps",sectionContent.str());
}

void THTMLPedestal::createNoiseDistribution()
{
	stringstream sectionContent;
	sectionContent<<"<h2> Noise Distribution of each detector</h2>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"Distribution of (ADC-Pedestal) for all events and all channels\n";
	sectionContent<< " This Distribution tells you what is the mean noise of each Detector\n";
//	sectionContent>>"";
	sectionContent<<"</p>\n";
	stringstream path;
	path<<this->path<<"/";
	sectionContent<<putImagesOfAllDetectors(path.str(),"hNosiseDistributionOfAllNonHitChannels_");
	this->addSection("Non-Hit Noise Distribution",sectionContent.str());
}

void THTMLPedestal::createHitOrderSection()
{
	stringstream sectionContent;
	sectionContent<<"<h2>Hit Order of each Detector Plane</h2>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"Order of highest signal to next adjacent highest signal. \n";
	sectionContent<< "this plot helps to figure out if there was any problem with \n";
	sectionContent<< "the Readout. A strong imbalance shows a dependency of readout direction.";
	sectionContent<<"</p>\n";
	stringstream path;
	path<<this->path<<"/";
	sectionContent<<putImagesOfAllDetectors(path.str(),"hSecondBiggestHitMinusBiggestHitPosition_");
	this->addSection("Hit Order",sectionContent.str());
}

void THTMLPedestal::createSaturatedChannels()
{

	stringstream sectionContent;
	sectionContent<<"<h2>Saturated Channels</h2>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"Histogramm of all Channels to see how often each channel got saturated\n";
//	sectionContent>>"";
	sectionContent<<"</p>\n";
	stringstream path;
	path<<this->path<<"/";
	sectionContent<<putImagesOfAllDetectors(path.str(),"hSaturatedChannels_");
	this->addSection("Saturated Channels",sectionContent.str());
}









