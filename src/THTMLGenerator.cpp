/*
 * THTMLGenerator.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: bachmair
 */

#include "../include/THTMLGenerator.hh"
using namespace std;
THTMLGenerator::THTMLGenerator(TSettings* newSettings) {
	// TODO Auto-generated constructor stub
	this->subdirPath="";
	this->mainPath="";
	this->path="";

	cout<<"GENERATE HTML FILE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
//	cout<<"GENERATE HTML FILE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
//	cout<<"GENERATE HTML FILE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
	settings=newSettings;
	if (settings==0)
		cerr<<"settings does not exist"<<endl;
	verbosity=3;
	sys= new TSystem();
	title = "Summary";
	setFileName("index.html");
//	cout<<"GENERATE HTML FILE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<fileName<<endl;
}

THTMLGenerator::~THTMLGenerator() {
	// TODO Auto-generated destructor stub
}


void THTMLGenerator::generateHTMLFile(){
	if (this->verbosity)cout<<"generateHTMLFile"<<endl;
	stringstream htmlOutputFileName;
	cout<<"FILE: "<<fileName<<endl;
	//cout<<"PATH: "<<sys->pwd()<<endl;
	htmlOutputFileName<<path<<"/"<<//sys->pwd()<<
			//"/"<<
			fileName;
	cout<<"create HTML file: \""<<htmlOutputFileName.str()<<"\""<<endl;
	html_summary.open(htmlOutputFileName.str().c_str());
	generatorHTMLHeader();
	generateHTMLTail();
	if (verbosity>2)cout<<"::GenerateHTML()::close html_summary"<<endl;
	html_summary.close();
}

void THTMLGenerator::generatorHTMLHeader()
{
		//summary page
	if (verbosity>2)cout <<"GENERATE HEADER"<<endl;
		html_summary << "<html>" << endl;
//		if (verbosity>2)cout << "Clustering::GenerateHTML():start html" <<eventReader<< endl;
		//Set Title

		html_summary << "<title>Run "<<settings->getRunNumber()<<" analysis results - "<<title<<"</title>" << endl;
		if (verbosity>2)cout << "Clustering::GenerateHTML():done with eventReader" << endl;

		html_summary << "<body bgcolor=\"White\">" << endl;
		html_summary << "<font face=\"Arial,Sans\" size=2>"<<endl;
		html_summary << "<a name=\"top\"></a>" << endl;

		if (verbosity>2)cout<<"Clustering::GenerateHTML():start summary"<<endl;
		html_summary << "<b>"<<title<<"</b>||"
				<< "<a href=\""<<this->mainPath<<"/index.html\">Summary</a>||"
				<< "<a href=\""<<this->mainPath<<"/pedestalAnalysis/pedestal.html\">Pedestal</a>||"
				<< "<a href=\""<<this->mainPath<<"/clustering/clustering.html\">Clustering</a>||"
				<< "<a href=\""<<this->mainPath<<"/selections/selection.html\">Selection</a>||"
				<< "<a href=\""<<this->mainPath<<"/selectionAnalysis/landaus.html\">Landaus</a>||"
//				<< "<a href=\"d8.html\">Diamond</a>||"
//				<< "<a href=\"d0.html\">D0X</a>||"
//				<< "<a href=\"d1.html\">D0Y</a>||"
//				<< "<a href=\"d2.html\">D1X</a>||"
//				<< "<a href=\"d3.html\">D1Y</a>||"
//				<< "<a href=\"d4.html\">D2X</a>||"
//				<< "<a href=\"d5.html\">D2Y</a>||"
//				<< "<a href=\"d6.html\">D3X</a>||"
//				<< "<a href=\"d7.html\">D3Y</a>"
				<<"<p>"
				<<endl;
		html_summary << "<center><font color=\"#000000\"><h1>Run "<<settings->getRunNumber()<<" analysis results - Summary</h1></font></center>"<<endl;
		html_summary << "<hr size=\"10\" Color=\"#ffff33\">" << endl;
		html_summary << "File generated at " << dateandtime.GetMonth() << "/ " << dateandtime.GetDay() << "/" << dateandtime.GetYear() << " at " << dateandtime.GetHour() << ":" << dateandtime.GetMinute() << ":" << dateandtime.GetSecond() << "<p>" << endl;
}

void THTMLGenerator::setFileName(string Name)
{
	fileName.clear();
	fileName="";
	stringstream output;
	output<<Name;
	size_t found= Name.find(".html");
	if(found==string::npos)
		output<<".html";
	fileName=output.str();
	cout<<"FileName set to \""<<fileName<<endl;
}

void THTMLGenerator::addSection(string sectionName, string secContent)
{
	cout<<"ADD SECTION NO"<<tableOfContent.size()<<":"<<endl;
	stringstream output;
	output<<"\t<h1 id=\"C"<<this->content.size()<<"\"><a name=\"C"<<this->content.size()<<"\">"<<sectionName<<"</a></h1>\n";
	output<<"\t<p>";
	output<<"\t"<<secContent;
	output<<"\t</p>";
//	cout<<output.str()<<endl;
	tableOfContent.push_back(sectionName);
	this->content.push_back(output.str());
}

void THTMLGenerator::generateTableOfContent(){
	if (verbosity>2)cout <<"GENERATE TABLE OF CONTENT"<<endl;

	stringstream output;
	output<<"\t<p>\n";
	for(UInt_t nSection=0;nSection<this->content.size();nSection++)
		output<<"\t<a href=\"#C"<<nSection<<"\">See also "<<tableOfContent.at(nSection)<<".</a><br />\n";
	output<<"\t</p>";
	if(verbosity>2)cout<<output.str()<<endl;
	html_summary<<output.str();
}
void THTMLGenerator::fillContent(){
	if (verbosity>2)cout <<"FILL CONTENT"<<endl;
	for(UInt_t nSec=0;nSec<content.size();nSec++){
		html_summary<<content.at(nSec);
	}
}
void THTMLGenerator::generateHTMLTail(){
	generateTableOfContent();
	fillContent();
	if (verbosity>2)cout <<"GENERATE HTML-TAIL"<<endl;
	html_summary << "<hr size=\"5\">" << endl;

	html_summary << "</font>"<<endl;
	html_summary<<"This Page was created on " << dateandtime.GetMonth() << "/ "
											  << dateandtime.GetDay() << "/"
											  << dateandtime.GetYear()
									<< " at " << dateandtime.GetHour() << ":"
										 	  << dateandtime.GetMinute() << ":"
										 	  << dateandtime.GetSecond()
				  << " with SVN Version No. " << SVN_REV<<"<p>";
	html_summary<<"</body>"<<endl;
	html_summary << "</HTML>" << endl;

}

void THTMLGenerator::setPathName(std::string pathName){
	this->path = pathName;
	cout<<"new Path: \""<<path<<"\""<<endl;
}

void THTMLGenerator::setMainPath(std::string mainPathName){
	this->mainPath=mainPathName;
	updatePath();
}

void THTMLGenerator::setSubdirPath(std::string subdirPathName){
	this->subdirPath=subdirPathName;
	updatePath();
}


std::string THTMLGenerator::createTable(std::vector<std::vector<std::string> > content)
{
	UInt_t nRows=content.size();
	UInt_t nCols =0;
	for(UInt_t row=0;row< nRows;row++)
		if(content.at(row).size()>nCols)nCols=content.at(row).size();
	cout<<"creating a Table with "<<nRows<<" Rows and "<<nCols<<" Columns!"<<endl;
	stringstream output;
	output<<"<p><table frame=\"void\" border=\"1\" rules=\"all\">\n";
	for(UInt_t row=0;row<nRows;row++){
		output<<"<tr>   ";
		for(UInt_t col=0;col<content.at(row).size();col++){
			//todo "Umlaute anpassen"
			if(row==0||col==0)	output<<"<th>"<<content.at(row).at(col)<<"</th>\t";
			else 		output<<"<td>"<<content.at(row).at(col)<<"</td>\t";
		}
		output<<"</tr>\n";
	}
	output<<"</table></p>\n";
	return (output.str());
}


std::string THTMLGenerator::putImage(std::string path, std::string name, std::string type,int percent)
{
	stringstream imageLink;
	imageLink<<path<<"//"<<name<<"."<<type;
	stringstream output;
	output<<"<img src=\""<<imageLink.str()<<"\" width=\""<<percent<<"%\" alt=\""<<name<<"\">\n";

	return (putLink(imageLink.str(),output.str()));
}
std::string THTMLGenerator::putLink(std::string link,std::string content){

	stringstream output;
	output<<"<a href=\""<<link<<"\">"<<content<<"</a>";
	return (output.str());
}

void THTMLGenerator::updatePath()
{
	setPathName(mainPath+subdirPath+"/");

}

std::string THTMLGenerator::floatToString(Float_t value)
{
	std::ostringstream out;
	out <<" "<< value<<" ";
	return (out.str());
}


std::string THTMLGenerator::putImagesOfAllDetectors(std::string path,std::string name, std::string type,int percentage){

	stringstream output;
	output<<"\n\t";
	for(UInt_t det = 0; det< TPlaneProperties::getNSiliconDetectors();det+=2){
		stringstream name2;
		name2<<name<<TADCEventReader::getStringForDetector(det);
		output<<putImage(path,name2.str());
	}
	output<<"\n<br\n\t";
	for(UInt_t det = 1; det< TPlaneProperties::getNSiliconDetectors();det+=2){
		stringstream name2;
		name2<<name<<TADCEventReader::getStringForDetector(det);
		output<<putImage(path,name2.str());
	}
	output<<"\n<br>\n\t";
	stringstream name2;
	name2<<name<<TADCEventReader::getStringForDetector(TPlaneProperties::getDetDiamond());
	output<<putImage(path,name2.str());
	return (output.str());
}
