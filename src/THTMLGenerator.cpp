/*
 * THTMLGenerator.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: bachmair
 */

#include "../include/THTMLGenerator.hh"
using namespace std;
THTMLGenerator::THTMLGenerator(TSettings* newSettings) {
	this->subdirPath="";
	this->mainPath="";
	this->path="";
	settings=newSettings;
	verbosity=settings->getVerbosity();
	if(verbosity)cout<<"\nGENERATE HTML FILE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;

	if (settings==0)
		cerr<<"settings does not exist"<<endl;
	setTitle("Summary");
	setFileName("overview.html");
}

THTMLGenerator::~THTMLGenerator() {
	if(verbosity)cout<<"delete HTML GEN"<<endl;
}



void THTMLGenerator::generateHTMLFile(){
	if (this->verbosity)cout<<"generateHTMLFile"<<endl;
	stringstream htmlOutputFileName;
	if(verbosity>2)cout<<"FILE: "<<fileName<<endl;
	htmlOutputFileName<<fileGenPath<<"/"<<	fileName;
	if(verbosity)cout<<"create HTML file: \""<<htmlOutputFileName.str()<<"\""<<endl;
	html_summary.open(htmlOutputFileName.str().c_str());
	generatorHTMLHeader();
	generateHTMLTail();
	if (verbosity>2)cout<<"::GenerateHTML()::close html_summary"<<flush;
	html_summary.close();
	if(verbosity>2)cout<<"...DONE"<<endl;
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
	string correctionPath;
	if(this->subdirPath!="")
		correctionPath="../";
	else correctionPath="./";

	html_summary << "<b>"<<title<<"</b>||"
			<< "<a href=\""<<correctionPath<<"overview.html\">Summary</a>||"
			<< "<a href=\""<<correctionPath<<"pedestalAnalysis/pedestal.html\">Pedestal</a>||"
			<< "<a href=\""<<correctionPath<<"clustering/clustering.html\">Clustering</a>||"
			<< "<a href=\""<<correctionPath<<"selections/selection.html\">Selection</a>||"
			<< "<a href=\""<<correctionPath<<"alignment/alignment.html\">Alignment</a>||"
			<< "<a href=\""<<correctionPath<<"selectionAnalysis/landaus.html\">Landaus</a>||"
            << "<a href=\""<<correctionPath<<"selectionAnalysis/selectionAnalysis.html\">selAnalysis</a>||"
			<< "<a href=\""<<correctionPath<<"transparentAnalysis/transparentAnalysis.html\">TransparentAnalysis</a>||"
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
	if(verbosity)cout<<"FileName set to \""<<fileName<<endl;
}

void THTMLGenerator::addSection(string sectionName, string secContent)
{
	if (verbosity>2)cout<<"ADD SECTION NO"<<tableOfContent.size()<<":"<<endl;
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
	if (verbosity>2)cout<<"new Path: \""<<path<<"\""<<endl;
}

void THTMLGenerator::setMainPath(std::string mainPathName){
	this->mainPath=mainPathName;
	updatePath();
}

void THTMLGenerator::setSubdirPath(std::string subdirPathName){
	this->subdirPath=subdirPathName;
	updatePath();
}

/**
 * create a string for a html table
 * @param content, outer vector: row,inner vector columm
 * @return html code for the table
 */
std::string THTMLGenerator::createTable(std::vector<std::vector<std::string> > content)
{
	UInt_t nRows=content.size();
	UInt_t nCols =0;
	for(UInt_t row=0;row< nRows;row++)
		if(content.at(row).size()>nCols)nCols=content.at(row).size();
	if (verbosity>2)cout<<"creating a Table with "<<nRows<<" Rows and "<<nCols<<" Columns!"<<endl;
	stringstream output;
	output<<"<p><table frame=\"void\" border=\"1\" rules=\"all\" align=\"center\">\n";
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


std::string THTMLGenerator::putImage(std::string path, TString name, std::string type,int percent,bool highlighted)
{
	std::string stringname;
	stringname.append(name);
	return putImage(path,stringname,type,percent);
}
std::string THTMLGenerator::putImage(std::string path, std::string name, std::string type,int percent,bool highlighted)
{
	stringstream imageLink;
	imageLink<<path<<"//"<<name<<"."<<type;
	stringstream output;
	output<<"<img src=\""<<imageLink.str()<<"\" width=\""<<percent<<"%\" alt=\""<<name<<"\"";
	if (highlighted)
	    output << " border=\"4\" ";
	output <<">\n";
	return (putLink(imageLink.str(),output.str()));
}
std::string THTMLGenerator::putLink(std::string link,std::string content){

	stringstream output;
	output<<"<a href=\""<<link<<"\">"<<content<<"</a>";
	return (output.str());
}

void THTMLGenerator::updatePath()
{
	if(verbosity)cout<<"THTMLGenerator::updatePath:";
	setPathName(mainPath+subdirPath+"/");
	if(verbosity)cout<<" \""<<this->path<<"\""<<endl;
}
string THTMLGenerator::putIntoCommand(string command, string input){
	stringstream output;
	output<<"<"<<command<<">"<<input<<"</"<<command<<">";
	return output.str();
}
string THTMLGenerator::center(string input)
{
	return putIntoCommand("center",input);
}

std::string THTMLGenerator::floatToString(Float_t value, UInt_t precision)
{
	std::ostringstream out;
	out <<" "<<std::setprecision(precision)<<value<<" ";
	return (out.str());
}


std::string THTMLGenerator::putImagesOfAllDetectors(std::string path,std::string name, std::string ending, std::string type,int percentage){
	cout<<"add all images for "<<name<<"DXX"<<ending<<" at path: "<<path<<endl;
	stringstream output;
	output<<"\n\t";
	for(UInt_t det = 0; det< TPlaneProperties::getNSiliconDetectors();det+=2){
		stringstream name2;
		name2<<name<<TPlaneProperties::getStringForDetector(det)<<ending;
		output<<putImage(path,name2.str());
	}
	output<<"\n<br>\n\t";
	for(UInt_t det = 1; det< TPlaneProperties::getNSiliconDetectors();det+=2){
		stringstream name2;
		name2<<name<<TPlaneProperties::getStringForDetector(det)<<ending;
		output<<putImage(path,name2.str());
	}
	output<<"\n<br>\n\t";
	stringstream name2;
	name2<<name<<TPlaneProperties::getStringForDetector(TPlaneProperties::getDetDiamond())<<ending;
	output<<putImage(path,name2.str());
	return (output.str());
}

string THTMLGenerator::cleanFileName(string fileName) {
//    while (fileName.find("//") !=fileName.end){
//        fileName.replace("//","");
//    }

    return fileName;
}
