/*
 * TRunInfo.cpp
 *
 *  Created on: Jul 17, 2012
 *      Author: bachmair
 */

#include "../include/TRunInfo.hh"

TRunInfo::TRunInfo() {
	// TODO Auto-generated constructor stub

}

TRunInfo::~TRunInfo() {
	// TODO Auto-generated destructor stub
}


std::string TRunInfo::getOutputDir()
{
	return outputDir;
}

void TRunInfo::setInputDir(std::string inputDir)
{
//	char resolved_path[200];
//	realpath(inputDir.c_str(), resolved_path);
	this->inputDir=getAbsPath(inputDir);
	cout<<"\nINPUTDIR: \""<<this->inputDir<<"\"\n";
}

void TRunInfo::setOutputDir(std::string outputDir)
{
//	char resolved_path[200];
//	realpath(outputDir.c_str(), resolved_path);

	this->outputDir = getAbsPath(outputDir);
	cout<<"\nOUTPUTDIR: \""<<this->outputDir<<"\"\n";
}

void TRunInfo::setRunSettingsDir(string settingsDir)
{
	this->runSettingsDir = getAbsPath(settingsDir);
	cout<<"\nsettingsDir: \""<<this->runSettingsDir<<"\"\n";
//	if(!this->runSettingsDir.compare("blaaaa"))
//		cout<<"failed..."<<endl;
//	else{
//		cout<<"yayyyy :D"<<endl;
//	}
//	char resolved_path[200];cout<<"\nbla0\n";
//	realpath(settingsDir.c_str(), resolved_path);cout<<"\nbla1\n";
//	printf("\nsettingsDir: \"%s\"\n",resolved_path);cout<<"\nbla2\n";
//	this->runSettingsDir = resolved_path;printf("\nblaaa3\n");
}

std::string TRunInfo::getAbsPath(std::string dir){
	char *resolved_path;
	resolved_path = realpath(dir.c_str(), NULL);
	if (resolved_path != NULL){
		string retrn_string(resolved_path);
		free(resolved_path);
		return retrn_string;
	}
	else{
		cout<<"It Failed..."<<endl;
		string bla("blaaaa");
		return bla;
	}
}





