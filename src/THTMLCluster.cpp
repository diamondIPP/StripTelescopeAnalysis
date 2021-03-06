/*
 * THTMLCluster.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: bachmair
 */

#include "../include/THTMLCluster.hh"

THTMLCluster::THTMLCluster(TSettings *settings):THTMLGenerator(settings) {
	setTitle("Clustering");
	this->setMainPath("../");
	this->setSubdirPath("clustering");
	this->setFileName("clustering.html");
	this->updatePath();
}

THTMLCluster::~THTMLCluster() {
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
		tablecontent.at(0).push_back(TPlaneProperties::getStringForDetector(det));
		tablecontent.at(1).push_back(floatToString(settings->getClusterSeedFactor(det,0)));
		tablecontent.at(2).push_back(floatToString(settings->getClusterHitFactor(det,0)));
	}
	sectionContent<<"<br>"<<this->createTable(tablecontent)<<"<br><br>";
	stringstream path;
	path<<this->path<<"/clustering/";
	for(UInt_t det = 0; det< TPlaneProperties::getNDetectors();det++){
		stringstream name;
		name<<"h2ndBiggestHitOverCharge_"<<TPlaneProperties::getStringForDetector(det);
		sectionContent<<putImage(path.str(),name.str());
	}
	this->addSection("Cluster Cuts",sectionContent.str());
}

void THTMLCluster::createContent()
{
	//	createClusterSize();
	createEtaDistributions();
}

void THTMLCluster::createEtaDistributions()
{
	stringstream sectionContent;

	//	sectionContent<<"<h2>Eta Distributions </h2>\n";
	sectionContent<<putImagesOfAllDetectors(".", "hEtaDistribution_");

	sectionContent<<"\n<br>\n\t";
	Int_t nDias = settings->getNDiamonds();
	for ( Int_t dia=0; dia<=nDias; dia++ ){
	    TString name =TString::Format("c_hEtaDistribution_Dia_Area%d",dia);
	    bool highlighted;
	    if (settings->getNDiamonds() ==2)
	    	highlighted = (dia ==  settings->getAnalysedDiamond());
	    else
	    	highlighted = (dia ==  settings->getAnalysedDiamond()+1);
	    sectionContent<<putImageOfPath((string)name,"png",20,highlighted );
	}
	sectionContent<<"\n<br>\n";
	sectionContent<<putImageOfPath("hEtaDistributionDia","png",20,false);

	this->addSection("Eta Distributions",sectionContent.str());
	sectionContent.clear();
	sectionContent.str("");
	sectionContent<<putImagesOfAllDetectors(".", "hEtaIntegral_");
	this->addSection("Eta Integrals",sectionContent.str());


	//c_hAsymmetricEtaFinal_
	sectionContent.clear();
	sectionContent.str("");
	sectionContent<<putImagesOfAllDetectors(".", "c_hAsymmetricEtaFinal_","All");
	sectionContent<<"\n<br>\n\t";
	for ( Int_t dia=0; dia<=nDias; dia++ ){
	    TString name =TString::Format("c_hAsymmetricEtaFinal_DiaArea%d",dia);
	    bool highlighted;
	    if (settings->getNDiamonds() ==2)
	    	highlighted = (dia ==  settings->getAnalysedDiamond());
	    else
	    	highlighted = (dia ==  settings->getAnalysedDiamond()+1);
	    sectionContent<<putImageOfPath((string)name,"png",20,highlighted );
	}
	sectionContent<<"\n<br>\n";
	sectionContent<<putImageOfPath("hAsymmetricEtaFinal_Dia","png",20,false);
	this->addSection("Cross Talk Corrected Eta Distributions",sectionContent.str());
}

void THTMLCluster::createPulseHeightPlots(vector<double> meanPulseHeigths)
{
	stringstream sectionContent;
	sectionContent<<"<p>\n"<<
			"Pulse Height Distribution for all ..."
			<<"</p>\n";
	std::vector< std::vector< std::string> > vecTable;
	if(meanPulseHeigths.size()<TPlaneProperties::getNDetectors()) meanPulseHeigths.resize(TPlaneProperties::getNDetectors());
	vecTable.resize(2);
	vecTable.at(0).push_back("");
	vecTable.at(1).push_back("mean PulseHeigth");
	for(UInt_t det=0;det <TPlaneProperties::getNDetectors();det++){
		vecTable.at(0).push_back(TPlaneProperties::getStringForDetector(det));
		vecTable.at(1).push_back(this->floatToString(meanPulseHeigths.at(det)));
	}
	sectionContent<<createTable(vecTable);
	sectionContent<<"\n\n<br><br>\n\n";
	sectionContent<<this->putImagesOfAllDetectors(".", "cPulseHeightDistribution_allClusterSizes_");
	addSection("Pulse Height Distribution",sectionContent.str());
}

void THTMLCluster::createClusterSize(std::vector<double> clusterSizes,std::vector<double> clusterSeedSizes,std::vector<double> numberOfClusters)
{
	stringstream sectionContent;
	std::vector< std::vector< std::string> > vecTable;
	vecTable.resize(4);
	vecTable.at(0).push_back("");
	vecTable.at(1).push_back("ClusterSize");
	vecTable.at(2).push_back("ClusterSeedSize");
	vecTable.at(3).push_back("Number of Clusters");
	if(clusterSizes.size()<TPlaneProperties::getNDetectors()) clusterSizes.resize(TPlaneProperties::getNDetectors());
	if(clusterSeedSizes.size()<TPlaneProperties::getNDetectors()) clusterSeedSizes.resize(TPlaneProperties::getNDetectors());
	if(numberOfClusters.size()<TPlaneProperties::getNDetectors()) numberOfClusters.resize(TPlaneProperties::getNDetectors());
	for(UInt_t det=0;det <TPlaneProperties::getNDetectors();det++){
		vecTable.at(0).push_back(TPlaneProperties::getStringForDetector(det));
		vecTable.at(1).push_back(this->floatToString(clusterSizes.at(det)));
		vecTable.at(2).push_back(this->floatToString(clusterSeedSizes.at(det)));
		vecTable.at(3).push_back(this->floatToString(numberOfClusters.at(det)));
	}
	sectionContent<<this->createTable(vecTable)<<"<br>\n";
	sectionContent<<"<h2>ClusterSize</h2>\n";
	sectionContent<<"<p>\n";
	//	sectionContent<<"Histogramm of all Channels to see how often each channel got saturated\n";
	//	sectionContent>>"";
	sectionContent<<"</p>\n";
	stringstream path;
	path<<this->path<<"/";
	stringstream fileNames;


	stringstream output,output2;
	output<<"\n\t";
	output2<<"\n\t";
	for(UInt_t det = 0; det< TPlaneProperties::getNSiliconDetectors();det+=2){
		Int_t seed = settings->getClusterSeedFactor(det,0);
		Int_t hit =settings->getClusterHitFactor(det,0);
		TString name =TString::Format("hClusterSize_Seed%d_Hit%d_",(int)seed,(int)hit).Append((string)TPlaneProperties::getStringForDetector(det));
//		name<<"hClusterSize_Seed"<<settings->getClusterSeedFactor(det,0)<<"_Hit"<<<<"_"<<;
		output<<putImage(path.str(),(string)name);
		name = TString::Format("hClusterSeedSize_Seed%d_Hit%d_",(int)seed,(int)hit).Append((string)TPlaneProperties::getStringForDetector(det));
		output2<<putImage(path.str(),(string)name);
	}
	output<<"\n<br\n\t";
	for(UInt_t det = 1; det< TPlaneProperties::getNSiliconDetectors();det+=2){
        Int_t seed = settings->getClusterSeedFactor(det,0);
        Int_t hit =settings->getClusterHitFactor(det,0);
        TString name =TString::Format("hClusterSize_Seed%d_Hit%d_",(int)seed,(int)hit).Append((TString)TPlaneProperties::getStringForDetector(det));
        output<<putImage(path.str(),(string)name);
        name = TString::Format("hClusterSeedSize_Seed%d_Hit%d_",(int)seed,(int)hit).Append((TString)TPlaneProperties::getStringForDetector(det));
        output2<<putImage(path.str(),(string)name);
	}
	output<<"\n<br>\n\t";
	output2<<"\n<br>\n\t";
	Int_t seed = settings->getClusterSeedFactor(TPlaneProperties::getDetDiamond(),0);
	Int_t hit = settings->getClusterHitFactor(TPlaneProperties::getDetDiamond(),0);
	TString name = TString::Format("hClusterSize_Seed%d_Hit%d_",seed,hit);
	name.Append((TString)TPlaneProperties::getStringForDetector(TPlaneProperties::getDetDiamond()));
	output<<putImage(path.str(),(string)name);

//	name<<"hClusterSize_Seed"<<seed<<"-Hit"<<hit<<"_"<<TPlaneProperties::getStringForDetector(TPlaneProperties::getDetDiamond());
//	name2<<"hClusterSeedSize_Seed"<<settings->getClusterSeedFactor(TPlaneProperties::getDetDiamond(),0)<<"-Hit"<<settings->getClusterHitFactor(TPlaneProperties::getDetDiamond(),0)<<"_"<<TPlaneProperties::getStringForDetector(TPlaneProperties::getDetDiamond());
	name = TString::Format("hClusterSeedSize_Seed%d_Hit%d_",seed,hit);
	name.Append((TString)TPlaneProperties::getStringForDetector(TPlaneProperties::getDetDiamond()));
	output2<<putImage(path.str(),(string)name);

	sectionContent<<output.str()<<"\n<br>\n";
	sectionContent<<output2.str()<<"\n<br>\n";

	sectionContent<<"<h2>NumberOfClusters </h2>\n";
	sectionContent<<putImagesOfAllDetectors(".", "NumberOfClusters_");
	this->addSection("Clusters",sectionContent.str());
}





