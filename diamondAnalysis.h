#include <fstream>
#include <iostream>
#include "TRunInfo.hh"

bool RunListOK;
std::string inputDir;
std::string outputDir;
std::string runSettingsDir="./";
std::string runListPath="RunList.ini";

int ReadRunList();

std::vector<TRunInfo> RunParameters;
