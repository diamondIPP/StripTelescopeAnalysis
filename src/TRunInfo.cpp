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
    this->inputDir = inputDir;
}

void TRunInfo::setOutputDir(std::string outputDir)
{
    this->outputDir = outputDir;
}



