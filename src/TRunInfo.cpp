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

std::string TRunInfo::getInputDir()
{
    return inputDir;
}

std::string TRunInfo::getOutputDir()
{
    return outputDir;
}

void TRunInfo::setInputDir(std::string inputDir)
{
  char resolved_path[300];
      realpath(inputDir.c_str(), resolved_path);
  printf("\nINPUTDIR: \"%s\"\t\n",resolved_path);
  this->inputDir=resolved_path;
}

void TRunInfo::setOutputDir(std::string outputDir)
{
  char resolved_path[300];
  realpath(outputDir.c_str(), resolved_path);
  printf("\nOUTPUTDIR: \"%s\"\n",resolved_path);

  this->outputDir = resolved_path;
}

void TRunInfo::setRunSettingsDir(string settingsDir)
{
  char resolved_path[300];
  realpath(settingsDir.c_str(), resolved_path);
  printf("\nsettingsDir: \"%s\"\n",resolved_path);
  this->runSettingsDir = resolved_path;
}





