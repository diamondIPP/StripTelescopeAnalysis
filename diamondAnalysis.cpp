/*
 * main.cpp
 * Diamond
 *
 * Created by Lukas Baeni on 19.01.11.
 */

//#include "TROOT.h"
//#include "Clustering.class.cpp"
#include "SlidingPedestal.class.hh"
#include "Clustering.class.hh"
#include "TRawEventSaver.hh"
#include "TPedestalCalculation.hh"
#include "TAnalysisOfClustering.hh"
#include "TClustering.hh"
#include "TSelectionClass.hh"
#include "diamondAnalysis.h"
#include "time.h"
#include "TSystem.h"
#include "TAlignment.hh"
#include "TSettings.class.hh"

using namespace std;
/*** USAGE ***/
void printHelp( void )
{
	cout<<"*******************************************************"<<endl;
	cout<<"diamondAnalysis: A Tool for anaylsis of RD42 test beams"<<endl;
	cout<<"USEAGE:\n"<<endl;
	cout<<"diamondAnalysis -i INPUTDIR -o OUTPUTDIR"<<endl;
	cout<<"*******************************************************"<<endl;

}
bool checkDir(string dir){
	return true;
}

bool readInputs(int argc,char ** argv){
	bool inputDirSet=false;
	bool outputDirSet=false;
	for(int i=1; i < argc; i++) {
			if(string(argv[i]) == "-h"||string(argv[i])=="--help")
					{
						printHelp();
						exit(0);
					}
		}
	for(int i=1; i<argc;i++) {
		if((string(argv[i]) == "-i"||string(argv[i])=="-I"||string(argv[i])=="--input"||string(argv[i])=="--INPUT")&&i+1<argc){
			i++;
			inputDir=string(argv[i]);
			inputDirSet=checkDir(inputDir);
			cout<<"found inputDir: \""<<inputDir<<"\""<<endl;
		}
		if((string(argv[i]) == "-o"||string(argv[i])=="-O"||string(argv[i])=="--output"||string(argv[i])=="--output")&&i+1<argc){
			i++;
			outputDir=string(argv[i]);
			outputDirSet=checkDir(outputDir);
			cout<<"found inputDir: \""<<outputDir<<"\""<<endl;
		}
		if((string(argv[i]) == "-r"||string(argv[i])=="-R")&&i+1<argc){
			i++;
			runListPath=string(argv[i]);
			cout<<"runListpath is set to:\""<<runListPath<<"\""<<endl;
		}

	}
	return true;
}



void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}


int main(int argc, char ** argv) {
	readInputs(argc,argv);
	cout<<"Currrent Subversion Revision: "<<SVN_REV<<endl;
	cout << "starting main loop.." << endl;
	initVariables();
	RunListOK = ReadRunList();
	TSystem* sys = gSystem;
	std::string currentDir = sys->pwd();
	cout << "Runnumbers ";
	for (unsigned int i = 0; i < RunParameters.size(); i++) {
		cout << RunParameters[i].RunNumber;
		if (i+1 < RunParameters.size()) cout << ", ";
	}
	cout << " will be analysed.." << endl;
	
	if (!RunListOK) return 0;
	
	/**Start with Analyising, read RunParameteres of the Run and start analysis with that parameters
	*/
	for (unsigned int i = 0; i < RunParameters.size(); i++) {
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
		cout << "ALTCLUSTERING: " << ALTCLUSTERING << endl;
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
		
		//FILE *log;
		
//		log = freopen(logfilename.str().c_str(), "w", stdout);

		//Save Events to RUNNUMBER/rawDATA.RUNNUMBER.root
		double vm2, rss2;
		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;


		stringstream settingsFileName;
		settingsFileName<<"settings."<<RUNNUMBER<<".ini";
		TSettings *settings=NULL;
		settings=new TSettings(settingsFileName.str(),RUNNUMBER);
		TRawEventSaver *eventSaver;
		eventSaver = new TRawEventSaver(RUNNUMBER);
		eventSaver->saveEvents(NEVENTS);
		delete eventSaver;

		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;


		//Calculate Pedestal
		sys->cd(currentDir.c_str());
		TPedestalCalculation* pedestalCalculation;
		pedestalCalculation = new TPedestalCalculation(RUNNUMBER,NEVENTS);
		pedestalCalculation->calculatePedestals(NEVENTS);
		pedestalCalculation->calculateSlidingPedestals(NEVENTS);
		delete pedestalCalculation;


		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;


		TClustering* clustering;
		clustering=new TClustering(RUNNUMBER);
		clustering->setSettings(settings);
		std::cout<<"cluster"<<endl;
		clustering->ClusterEvents(NEVENTS);
		delete clustering;

		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;

		TSelectionClass* selectionClass;
		selectionClass=new TSelectionClass(settings);
		selectionClass->MakeSelection(NEVENTS);
		delete selectionClass;


		TAlignment *alignment;
		alignment= new TAlignment(settings);
		alignment->setSettings(settings);
//		alignment->Align();
		delete alignment;

//
//		TAnalysisOfClustering* analysisClustering;
//		analysisClustering= new TAnalysisOfClustering(RUNNUMBER);
//		analysisClustering->doAnalysis(NEVENTS);
//		delete analysisClustering;

		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;


		if (settings!=NULL)delete settings;
//		if (DO_SLIDINGPEDESTAL) {
//			cout << endl;
//			cout << "==> Starting SlidingPedestal.." << endl;
//			cout << "SlidingPedestal sl(" << RUNNUMBER << ",\"" << RUNDESCRIPTION << "\");" << endl;
//			SlidingPedestal sl(RUNNUMBER,RUNDESCRIPTION);
//			cout << "sl.Slide(" << NEVENTS << "," << INITIAL_EVENT << "," << HIT_OCCUPANCY << ");" << endl;
//			sl.Slide(NEVENTS,INITIAL_EVENT,HIT_OCCUPANCY);
//		}
//		cout << endl;
//		cout << "==> Starting Clustering.." << endl;
//		cout << "Clustering cl(" << RUNNUMBER << ",\"" << RUNDESCRIPTION << "\")" << endl;
//		Clustering cl(RUNNUMBER,RUNDESCRIPTION);
//		cl.verbosity = VERBOSITY;
//		vector<FidCutRegion> FidCutRegions;
//		if (cl.getUseAutoFidCut()) {
//			cout << endl;
//			cout << "==> Starting AutoFidCut.." << endl;
//			cout << "cl.AutoFidCut()" << endl;
//			cl.AutoFidCut();
//			if (FidCutRegions.size() == 0) cl.setUseAutoFidCut(false);
//		}
//		if (FidCutRegions.size() > 0 && cl.getUseAutoFidCut()) {
//			for (int reg = 0; reg < FidCutRegions.size(); reg++) {
//				cl.SetRunParameters(reg,FidCutRegions[reg],FidCutRegions.size()-1);
//				// TODO: set different paths for the plots
//				cout << "cl.ClusterRun(" << PLOTS << ");" << endl;
//				cl.ClusterRun(PLOTS);
//			}
//		}
//		else {
//			cl.setAlternativeClustering(ALTCLUSTERING);
//			if (DO_ALIGNMENT) {
//				cout << "cl.Align(" << PLOTS << "," << CUTFAKETRACKS << ");" << endl;
//				cl.Alignment(PLOTS, CUTFAKETRACKS);
//			}
//			else {
//				cout << "cl.ClusterRun(" << PLOTS << ");" << endl;
//				cl.ClusterRun(PLOTS);
//			}
//		}
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
	ALTCLUSTERING = 0;
    VERBOSITY=0;

}

int ReadRunList() {
	RunInfo run;
	char RunDescription[200];
	int NEvents, Initial_Event;
	RunParameters.clear();
	cout << endl << "reading runlist.." << endl;
	ifstream file(runListPath.c_str());//"RunList.ini");
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
		
		
		sscanf(line.c_str(), "%d %s %d %d %d %d %d %d", &RUNNUMBER, RunDescription, &VERBOSITY, &NEvents, &Initial_Event, &CUTFAKETRACKS, &DO_SLIDINGPEDESTAL, &DO_ALIGNMENT/*, &ALTCLUSTERING*/);
		if (NEvents != 0) NEVENTS = NEvents;
		if (Initial_Event != 0) INITIAL_EVENT = Initial_Event;
		cout << "RunDescription Char: " << RunDescription[0] << endl;
		if (RunDescription[0] != '0') RUNDESCRIPTION = RunDescription;
		run.SetParameters();
		RunParameters.push_back(run);
	}
	return 1;
}
