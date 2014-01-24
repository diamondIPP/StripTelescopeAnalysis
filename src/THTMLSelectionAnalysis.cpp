/**
 * @file THTMLSelectionAnalysis.cpp
 *
 * @date Feb 12, 2013
 * @author bachmair
 * @description
 */

#include "../include/THTMLSelectionAnalysis.hh"

THTMLSelectionAnalysis::THTMLSelectionAnalysis(TSettings *settings):THTMLGenerator(settings)  {
	// TODO Auto-generated constructor stub
	this->setFileName("selectionAnalysis.html");
	this->setMainPath("../");
	this->setSubdirPath("selectionAnalysis/");
	this->setTitle("Selection Analysis");
}

THTMLSelectionAnalysis::~THTMLSelectionAnalysis() {
	// TODO Auto-generated destructor stub
}

void THTMLSelectionAnalysis::addSelectionPlots() {
	stringstream sectionContent;
	sectionContent<<"<p>\n";
	sectionContent<<this->putImage(path,(string)"chValidSiliconAndDiamondHit",(string)"png",32)<<"\n";
	sectionContent<<this->putImage(path,(string)"chValidSiliconAndOneDiamondHit",(string)"png",32)<<"\n";
	sectionContent<<this->putImage(path,(string)"chValidSiliconAndOneDiamondHitNotMaskedAdjacentChannels",(string)"png",32)<<"<br>\n";
	sectionContent<<this->putImage((string)path,(string)"chFidCut_oneDiamondCluster",(string)"png",49)<<"\n";

	sectionContent<<this->putImage(path,(string)"chValidSiliconAndOneDiamondHitInSameAreaAndFidCut",(string)"png",49)<<"<br>\n";
	sectionContent<<"</p>";
	this->addSection("Selection Plots",sectionContent.str());
}

void THTMLSelectionAnalysis::addAreaPlots() {
	stringstream sectionContent;
	sectionContent<<"<p>\n";
	for(Int_t area=0;area<settings->getNDiaDetectorAreas();area++){
		int firstCh = settings->getDiaDetectorArea(area).first;
		int lastCh = settings->getDiaDetectorArea(area).second;
		TString name = TString::Format("hChargeOfCluster_ClusterSize_1_2_2D_noBorderHit_area_%d_ch_%d_%d",area,firstCh,lastCh);
		sectionContent<<"<h2>"<<TString::Format("Area %d: Channel %d - %d ",area,firstCh,lastCh)<<"</h2>\n";
		sectionContent<<this->putImage(path,name,(string)"png",49)<<"\n";

		name = TString::Format("cChargeOfCluster_ClusterSize_1_2_NoBorderHit_area_%d_ch_%d_%d",area,firstCh,lastCh);
		sectionContent<<this->putImage(path,name,(string)"png",49)<<"<br>\n";
		name = TString::Format("cChargePerChannel_area%d",area);
		sectionContent<<this->putImage(path,name,(string)"png",33)<<"\n";
		name = TString::Format("hClusterSizeVsChannelPos_Area_%d_ch_%d_%d",area,firstCh,lastCh);
		sectionContent<<this->putImage(path,name,(string)"png",33)<<"\n";
		name = TString::Format("hClusterSize_Area_%d_ch_%d_%d",area,firstCh,lastCh);
		sectionContent<<this->putImage(path,name,(string)"png",33)<<"<br>\n";
		sectionContent<<"<br>\n";
	}
	sectionContent<<"</p>";

	this->addSection("Area Plots",sectionContent.str());
}

void THTMLSelectionAnalysis::addFiducialCutPlots() {
	stringstream sectionContent;
	sectionContent<<"<p>\n";
	for(UInt_t i=0;i<settings->getSelectionFidCuts()->getNFidCuts();i++){
		sectionContent<<"<h2>"<<TString::Format("Fiducial Cut %d:",i+1)<<"</h2>\n";

		TString name = TString::Format("hMeanChargeFiducialCutNo%d",i+1);
		sectionContent<<this->putImage(path,name,(string)"png",33)<<"\n";
		name = TString::Format("hChargeVsFidX_HitInFidCutNo%d",i+1);
		sectionContent<<this->putImage(path,name,(string)"png",33)<<"\n";
		name = TString::Format("hChargeVsFidY_HitInFidCutNo%d",i+1);
		sectionContent<<this->putImage(path,name,(string)"png",33)<<"\n";
		sectionContent<<"<br><br>\n";
	}
	sectionContent<<"</p>";
	this->addSection("FiducialCut Plots",sectionContent.str());
}
