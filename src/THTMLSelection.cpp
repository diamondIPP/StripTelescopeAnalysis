/*
 * THTMLSelection.cpp
 *
 *  Created on: May 16, 2012
 *      Author: bachmair
 */

#include "../include/THTMLSelection.hh"

THTMLSelection::THTMLSelection(TSettings *settings):THTMLGenerator(settings) {
	this->setFileName("selection.html");
	this->setSubdirPath("selections/");
	this->setMainPath("../");
	this->setTitle("Selection - Cut Flow");

	this->updatePath();
}

THTMLSelection::~THTMLSelection() {
	// TODO Auto-generated destructor stub
	cout<<"Delete HTMLSelection"<<endl;
}

void THTMLSelection::createCutFlowTable(std::vector<int> vecCutFlow	)
{
}


void THTMLSelection::createFiducialCuts()
{
	stringstream sectionContent;
	sectionContent<<"<h3>Fiducial Cut - valid Silicon Track</h3>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"The Fiducial cut is applied to the mean of the position of the cluster in each plane.\n ";
	sectionContent<<"To find good values for the fiducial cut. The mean position of all silicon planes is calculated\n ";
	sectionContent<<"and plotted in the next plot. To be able to calculate this mean position one and only one Cluster\n";
	sectionContent<<"in each detector is required.<br>\n";
	sectionContent<<(putImage(this->path,(string)"chFidCutSilicon_OneAndOnlyOneCluster","png",50));
	sectionContent<<"</p>\n";
	sectionContent<<"<br>\n";
	sectionContent<<"<h3>Fiducial Cut - valid Silicon Track && one Diamond Cluster</h3>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"The next plot shows the mean position of the clusters in each silicon plane with one and only one\n";
	sectionContent<<"cluster in the diamond detector. With this condition you can see the needed fiducial cuts to calculate\n";
	sectionContent<<"efficency in the diamond detector.<br>\n";
	sectionContent<<(putImage(this->path,(string)"chFidCutSilicon_OneAndOnlyOneCluster_DiamondCluster","png",50))<<"<br>\n";

	sectionContent<<putImage(this->path,(string)"chProjX","png",49)<<" \n";
	sectionContent<<putImage(this->path,(string)"chProjY","png",49)<<"<br>\n";
	sectionContent<<putImage(this->path,(string)"chSelectedEvents","png",50)<<" \n";

	sectionContent<<"</p>\n";
	this->addSection("Fiducial Cut",sectionContent.str());
}

void THTMLSelection::createCutFlowGraph(std::string content)
{
	size_t found = content.find_first_of("\n");
	while(found!=string::npos){
		content.replace(found,1,"<br>");
		found = content.find_first_of("\n");
	}
	stringstream sectionContent;
	sectionContent<<"<p>\n<pre>\n"<<content<<"\n</pre>\n</p>\n";

	sectionContent<<"<h3>Cut Flow Diagrams</h3>\n";
	sectionContent<<"<h4>Silicon Cuts</h4>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"First Pie Chart shows the cuts for the silicon Plane, <br>\n";
	sectionContent<<" - noSilTrack (red) means that there is not in each Silicon Detector one and only one Silicon Cluster<br>\n";
	sectionContent<<" - notInFidCut (green) means that there is one and only one Cluster in each silicon plane, but the mean Position is not lying in the fiducial cuts<br>\n";
	sectionContent<<" - notExactlyOneDiamondCluster (blue), means that it has a valid silicon Track but there is not exactly one Cluster in the diamond detector<br>\n";
	sectionContent<<" - useForAlignment (yellow) are Events which have a valid Sil Track and have one and only one Cluster in the diamond hit. With the training track fraction they are choosen for Alignment<br>\n";
	sectionContent<<" - useForAnalysis (pink) are Events which have a valid Sil Track and have one and only one Cluster in the diamond hit. With the training track fraction they are choosen for Analysis<br>\n";
	sectionContent<<putImage(this->path,(string)"cMainCutFlow","png",50)<<"<br>\n";
	sectionContent<<"</p>\n";
	sectionContent<<"<h4>Diamond Cuts</h4>\n";
	sectionContent<<"<p>\n";
	sectionContent<<"This is a Pie Chart of all Events which have a valid Silicon Track (e.g. one and only one Cluster in each detector && mean is in fiducialCut area).<br>\n";
	sectionContent<<"the point notExactlyOneDiamondCluster from the chart above is divided into moreThanOneDiamondCluster and noDiamondCluster<br>\n";
	sectionContent<<"with that pie chart you can see the efficency of the diamond <br>\n";
	sectionContent<<putImage(this->path,(string)"cPieValidSiliconTrack","png",50);
	sectionContent<<"<br>\n";

	sectionContent<<putImage(this->path,(string)"hAnalysisFraction","png",70);
	sectionContent<<"<br></p>";
	this->addSection("CutFlow",sectionContent.str());

}




