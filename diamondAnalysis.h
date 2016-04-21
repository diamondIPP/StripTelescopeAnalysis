#include <fstream>
#include <iostream>
#include "TRunInfo.hh"

bool RunListOK;
std::string inputDir="./";
std::string outputDir="./";
std::string runSettingsDir="./";
std::string runListPath="RunList.ini";
bool run_3danalysis = true;

int ReadRunList();

std::vector<TRunInfo> RunParameters;
