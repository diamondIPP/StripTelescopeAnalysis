/*
 * THTMLGenerator.hh
 *
 *  Created on: Feb 15, 2012
 *      Author: bachmair
 */

#ifndef THTMLGENERATOR_HH_
#define THTMLGENERATOR_HH_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TSystem.h"
#include "TDatime.h"
#include "TSettings.class.hh"
#include "TADCEventReader.hh"
#include "TDetectorAlignment.hh"


class THTMLGenerator {
public:
	THTMLGenerator(TSettings *settings);
	virtual ~THTMLGenerator();
	void generateHTMLFile();
	void setFileName(string Name);
	void setPathName(string pathName);
	void setMainPath(std::string mainPathName);
	void setFileGeneratingPath(string newFileGenPath){if(verbosity)cout<<"FileGeneration at:"<<newFileGenPath<<endl;fileGenPath=newFileGenPath;};
	void setSubdirPath(std::string subDirPath);
	void updatePath();
	void addSection(string sectionName, string content);
	void setTitle(string newTitle){if(verbosity)cout<<"HTML:SetTitle:"<<newTitle<<endl;this->title.assign(newTitle);};
	virtual void createContent(){};
public:
//	template <typename T1, typename T2>
//	string combineToString(T1 a, T2 b){
//	    stringstream output;
//	    output<<a<<b;
//	    return output.str();
//	  }
string combineToString(string a,UInt_t b){
  stringstream output;
  output<<a<<b;
  return output.str();
}
string combineToString(string a,Double_t b,UInt_t precision=3){
  stringstream output;
  output<<a<<std::setprecision(precision)<<b;
  return output.str();
}
protected:
 string putIntoCommand(string command, string input);
  string center(string input);
	void generatorHTMLHeader();
	void generateHTMLTail();
	void generateTableOfContent();
	void fillContent();
	std::string createTable(std::vector<std::vector<std::string> > content);
	std::string putImage(std::string path, std::string name, std::string type = "png",int sizeInPercentage=20);
public:
	std::string putImageOfPath(std::string name, std::string type="png",int sizeInPercentage=20){return putImage(path,name,type,sizeInPercentage);};
protected:
	std::string fileGenPath;
	std::string putImagesOfAllDetectors(std::string path,std::string name, std::string type="png", int percentage =20);
	std::string putLink(std::string link,std::string content);
	TSettings* settings;
	UInt_t verbosity;
	TSystem *sys;
	std::ofstream html_summary;
    TDatime dateandtime;
    std::string fileName;
    std::string path;
    std::string mainPath;
    std::string subdirPath;
    std::vector<std::string> content;
    std::vector<std::string> tableOfContent;
    std::string title;
public:
    static std::string floatToString(Float_t value, UInt_t precision=4);
};

#endif /* THTMLGENERATOR_HH_ */
