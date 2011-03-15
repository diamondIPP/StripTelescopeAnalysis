/*
 * main.cpp
 * Diamond
 *
 * Created by Lukas Baeni on 19.01.11.
 */

//#include "TROOT.h"
//#include "Clustering.class.cpp"
#include "SlidingPedestal.class.cpp"
#include "Clustering.class.cpp"
#include <fstream>
#include <iostream>
#include "main.h"
#include "time.h"

using namespace std;

int main () {
	cout << "starting main loop.." << endl;
	initVariables();
	RunListOK = ReadRunList();
	
	cout << "Runnumbers ";
	for (int i = 0; i < RunParameters.size(); i++) {
		cout << RunParameters[i].RunNumber;
		if (i+1 < RunParameters.size()) cout << ", ";
	}
	cout << " will be analysed.." << endl;
	
	if (!RunListOK) return 0;
	
	for (int i = 0; i < RunParameters.size(); i++) {
		RunParameters[i].GetParameters();
		cout << endl << endl << endl << endl;
		cout << "====================================" << endl;
		cout << "==> Starting analysis.." << endl;
		cout << "====================================" << endl << endl;
		cout << "RUNNUMBER: " << RUNNUMBER << endl;
		cout << "NEVENTS: " << NEVENTS << endl;
		cout << "RUNDESCRIPTION: " << RUNDESCRIPTION << endl;
		cout << "VERBOSITY: " << VERBOSITY << endl;
		cout << "INITIAL_EVENT: " << INITIAL_EVENT << endl;
		cout << "HIT_OCCUPANCY: " << HIT_OCCUPANCY << endl;
		cout << "ALTERNATIVECLUSTERING: " << ALTERNATIVECLUSTERING << endl;
		cout << "DO_ALIGNMENT: " << DO_ALIGNMENT << endl;
		cout << "DO_SLIDINGPEDESTAL: " << DO_SLIDINGPEDESTAL << endl;
		cout << "CUTFAKETRACKS: " << CUTFAKETRACKS << endl;
		cout << endl << endl << endl;
		
		time_t rawtime;
		tm *timestamp;
		
		time (&rawtime);
		
		timestamp = gmtime(&rawtime);
		
		ostringstream logfilename;
		logfilename << "analyse_log_" << RUNNUMBER << "_" << timestamp->tm_year << "-" << timestamp->tm_mon << "-" << timestamp->tm_mday << "." << timestamp->tm_hour << "." << timestamp->tm_min << "." << timestamp->tm_sec << ".log";
		
		FILE *log;
		
//		log = freopen(logfilename.str().c_str(), "w", stdout);
		
		if (DO_SLIDINGPEDESTAL) {
			cout << endl;
			cout << "==> Starting SlidingPedestal.." << endl;
			cout << "SlidingPedestal sl(" << RUNNUMBER << ",\"" << RUNDESCRIPTION << "\");" << endl;
			SlidingPedestal sl(RUNNUMBER,RUNDESCRIPTION);
			cout << "sl.Slide(" << NEVENTS << "," << INITIAL_EVENT << "," << HIT_OCCUPANCY << ");" << endl;
			sl.Slide(NEVENTS,INITIAL_EVENT,HIT_OCCUPANCY);
		}
		cout << endl;
		cout << "==> Starting Clustering.." << endl;
		cout << "Clustering cl(" << RUNNUMBER << ",\"" << RUNDESCRIPTION << "\")" << endl;
		Clustering cl(RUNNUMBER,RUNDESCRIPTION);
		cl.Verbosity = VERBOSITY;
		vector<FidCutRegion> FidCutRegions;
		if (cl.UseAutoFidCut) {
			cout << endl;
			cout << "==> Starting AutoFidCut.." << endl;
			cout << "cl.AutoFidCut()" << endl;
			cl.AutoFidCut();
			if (FidCutRegions.size() == 0) cl.UseAutoFidCut = false;
		}
		if (FidCutRegions.size() > 0 && cl.UseAutoFidCut) {
			for (int reg = 0; reg < FidCutRegions.size(); reg++) {
				cl.SetRunParameters(reg,FidCutRegions[reg],FidCutRegions.size()-1);
				// TODO: set different paths for the plots
				cout << "cl.ClusterRun(" << PLOTS << ");" << endl;
				cl.ClusterRun(PLOTS);
			}
		}
		else {
			cl.AlternativeClustering = ALTERNATIVECLUSTERING;
			if (DO_ALIGNMENT) {
				cout << "cl.Align(" << PLOTS << "," << CUTFAKETRACKS << ");" << endl;
				cl.Align(PLOTS, CUTFAKETRACKS);
			}
			else {
				cout << "cl.ClusterRun(" << PLOTS << ");" << endl;
				cl.ClusterRun(PLOTS);
			}
		}
//	    fclose(log);	
	}

	return 0;
}

void initVariables() {
	RUNDESCRIPTION = "";
	NEVENTS = 10000;
	INITIAL_EVENT = 1000;
	HIT_OCCUPANCY = 0;
	PLOTS = 1;
	ALTERNATIVECLUSTERING = 0;
    VERBOSITY=0;
}

int ReadRunList() {
	RunInfo run;
	char RunDescription[200];
	int NEvents, Initial_Event;
	RunParameters.clear();
	cout << endl << "reading runlist.." << endl;
	ifstream file("RunList.ini");
	if (!file) {
		cout << "An error has encountered while trying to open RunList.ini" << endl;
		return 0;
	}
	else cout << "RunList.ini" << " successfully opened." << endl << endl;
	
	while (!file.eof()) {
		initVariables();
		
		//get next line
		string line;
		getline(file,line);
		
		//check if comment or empty line
		if ((line.substr(0, 1) == ";") || (line.substr(0, 1) == "#") || (line.substr(0, 1) == "/") || line.empty()) {
			continue;
		}
		
		
		sscanf(line.c_str(), "%d %s %d %d %d %d %d %d %d", &RUNNUMBER, RunDescription, &VERBOSITY, &NEvents, &Initial_Event, &CUTFAKETRACKS, &DO_SLIDINGPEDESTAL, &DO_ALIGNMENT, &ALTERNATIVECLUSTERING);
		if (NEvents != 0) NEVENTS = NEvents;
		if (Initial_Event != 0) INITIAL_EVENT = Initial_Event;
		cout << "RunDescription Char: " << RunDescription[0] << endl;
		if (RunDescription[0] != '0') RUNDESCRIPTION = RunDescription;
		cout << "RUNNUMBER: " << RUNNUMBER << endl;
		cout << "NEVENTS: " << NEVENTS << endl;
		cout << "RUNDESCRIPTION: " << RUNDESCRIPTION << endl;
		cout << "VERBOSITY: " << VERBOSITY << endl;
		cout << "INITIAL_EVENT: " << INITIAL_EVENT << endl;
		cout << "HIT_OCCUPANCY: " << HIT_OCCUPANCY << endl;
		cout << "ALTERNATIVECLUSTERING: " << ALTERNATIVECLUSTERING << endl;
		cout << "DO_ALIGNMENT: " << DO_ALIGNMENT << endl;
		cout << "DO_SLIDINGPEDESTAL: " << DO_SLIDINGPEDESTAL << endl;
		cout << "CUTFAKETRACKS: " << CUTFAKETRACKS << endl;
		run.SetParameters();
		RunParameters.push_back(run);
	}
}
