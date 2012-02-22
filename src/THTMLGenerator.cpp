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
	settings=newSettings;
	if (settings==0)
		cerr<<"settings does not exist"<<endl;
	verbosity=1;
	sys= new TSystem();

}

THTMLGenerator::~THTMLGenerator() {
	// TODO Auto-generated destructor stub
}

void THTMLGenerator::createAlignmentSummary(TDetectorAlignment *alignment)
{

}

void THTMLGenerator::generateHTMLFile(){
	if (this->verbosity)cout<<"Clustering::GenerateHTML()"<<endl;
	stringstream htmlOutputFileName;
	htmlOutputFileName<<sys->pwd()<<"/index.html";
	cout<<"create HTML file: \""<<htmlOutputFileName.str()<<"\""<<endl;
	html_summary.open(htmlOutputFileName.str().c_str());
	generatorHTMLHeader();
}

void THTMLGenerator::generatorHTMLHeader()
{
		int section, subsection;

		//summary page

		html_summary << "<html>" << endl;
//		if (verbosity>2)cout << "Clustering::GenerateHTML():start html" <<eventReader<< endl;

		html_summary << "<title>Run "<<settings->getRunNumber()<<" analysis results - Summary</title>" << endl;
		if (verbosity>2)cout << "Clustering::GenerateHTML():done with eventReader" << endl;
		html_summary << "<body bgcolor=\"White\">" << endl;
		html_summary << "<font face=\"Arial,Sans\" size=2>"<<endl;
		html_summary << "<a name=\"top\"></a>" << endl;

		if (verbosity>2)cout<<"Clustering::GenerateHTML():start summary"<<endl;
		html_summary << "<b>Summary</b>||"
				<< "<a href=\"d8.html\">Diamond</a>||"
				<< "<a href=\"d0.html\">D0X</a>||"
				<< "<a href=\"d1.html\">D0Y</a>||"
				<< "<a href=\"d2.html\">D1X</a>||"
				<< "<a href=\"d3.html\">D1Y</a>||"
				<< "<a href=\"d4.html\">D2X</a>||"
				<< "<a href=\"d5.html\">D2Y</a>||"
				<< "<a href=\"d6.html\">D3X</a>||"
				<< "<a href=\"d7.html\">D3Y</a><p>"<<endl;
		html_summary << "<center><font color=\"#000000\"><h1>Run "<<settings->getRunNumber()<<" analysis results - Summary</h1></font></center>"<<endl;
		html_summary << "<hr size=\"10\" Color=\"#ffff33\">" << endl;
		html_summary << "Results from " << dateandtime.GetMonth() << "/ " << dateandtime.GetDay() << "/" << dateandtime.GetYear() << " at " << dateandtime.GetHour() << ":" << dateandtime.GetMinute() << ":" << dateandtime.GetSecond() << ")<p>" << endl;


}

void THTMLGenerator::generateHTMLTail(){

	html_summary << "<hr size=\"5\">" << endl;

	html_summary << "</font>"<<endl;
	html_summary << "</HTML>" << endl;
	if (verbosity>2)cout<<"Clustering::GenerateHTML()::close html_summary"<<endl;
	html_summary.close();
}


