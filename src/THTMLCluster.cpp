/*
 * THTMLCluster.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: bachmair
 */

#include "../include/THTMLCluster.hh"

THTMLCluster::THTMLCluster(TSettings *settings):THTMLGenerator(settings) {
	// TODO Auto-generated constructor stub

	setTitle("Clustering");

}

THTMLCluster::~THTMLCluster() {
	// TODO Auto-generated destructor stub
}

void THTMLCluster::createTableOfCuts()
{
	stringstream sectionContent;
	sectionContent<<"<h1>Seed and Hit Values in Units of Sigma</h1>\n";
	std::vector<std::vector< std::string > > tablecontent;
	tablecontent.resize(3);
	tablecontent.at(0).push_back("Detector");
	tablecontent.at(1).push_back("Seed");
	tablecontent.at(2).push_back("Hit");
	for(UInt_t det =0;det <TPlaneProperties::getNDetectors();det++){
		tablecontent.at(0).push_back(TADCEventReader::getStringForDetector(det));
		tablecontent.at(1).push_back(floatToString(settings->getClusterSeedFactor(det)));
		tablecontent.at(2).push_back(floatToString(settings->getClusterHitFactor(det)));
	}
	sectionContent<<"<br>"<<this->createTable(tablecontent)<<"<br><br>";
	stringstream path;
	path<<this->path<<"/clustering/";
	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
		stringstream name;
		name<<"h2ndBiggestHitOverCharge_"<<TADCEventReader::getStringForDetector(det);
		sectionContent<<putImage(path.str(),name.str());
	}
	this->addSection("Cluster Cuts",sectionContent.str());
}



