/*
 * main.cpp
 * Diamond
 *
 * Created by Lukas Baeni on 19.01.11.
 */

//#include "TROOT.h"
//#include "Clustering.class.cpp"
//#include "SlidingPedestal.class.hh"
//#include "Clustering.class.hh"
#include "TRawEventSaver.hh"
#include "TPedestalCalculation.hh"
#include "TAnalysisOfPedestal.hh"
#include "TAnalysisOfClustering.hh"
#include "TAnalysisOfSelection.hh"
#include "TClustering.hh"
#include "TSelectionClass.hh"
#include "THTMLGenerator.hh"
#include "THTMLCluster.hh"
#include "diamondAnalysis.h"
#include "time.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TAlignment.hh"
#include "TSettings.class.hh"
#include "TTransparentAnalysis.hh"
#include "TAnalysisOfAlignment.hh"
#include "TResults.hh"

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

  TStopwatch comulativeWatch;
  comulativeWatch.Start(true);
	readInputs(argc,argv);
	cout<<"Currrent Subversion Revision: "<<SVN_REV<<endl;
	cout << "starting main loop.." << endl;
	RunListOK = ReadRunList();
	TSystem* sys = gSystem;
//  if(gStyle!=0)
//    if(!gStyle->IsZombie()){
//      gROOT->SetStyle("Plain"); //General style (see TStyle)
////      gStyle->SetOptStat(221111111); //Stat options to be displayed     without under- and overflow use gStyle->SetOptStat(1110);
//      if(gStyle->GetOptStat()!=221111111)
//        gStyle->SetOptStat("nemrKSiou");
//      gStyle->SetOptFit(1111);  //Fit options to be displayed
//      gStyle->SetPadBottomMargin(0.15); //Gives more space between histogram and edge of plot
//      gStyle->SetPadRightMargin(0.15);
//      gStyle->SetPadTopMargin(0.15);
//      //gStyle->SetTitleColor(19,"");
//      gStyle->SetStatH(0.12); //Sets Height of Stats Box
//      gStyle->SetStatW(0.15); //Sets Width of Stats Box
//      gStyle->SetPalette(1); // determines the colors of temperature plots (use 1 for standard rainbow; 8 for greyscale)
//    }
	std::string currentDir = sys->pwd();
	for (unsigned int i = 0; i < RunParameters.size(); i++) {
		cout << RunParameters[i].getRunNumber();
		if (i+1 < RunParameters.size()) cout << ", ";
	}
	cout << " will be analysed.." << endl;
	
	if (!RunListOK) return 0;
	
	/**Start with Analyising, read RunParameteres of the Run and start analysis with that parameters
	*/
	for (unsigned int i = 0; i < RunParameters.size(); i++) {

	  TStopwatch runWatch;
	  runWatch.Start(true);
		UInt_t RUNNUMBER = RunParameters[i].getRunNumber();
		UInt_t VERBOSITY = RunParameters[i].getVerbosity();
		std::string RUNDESCRIPTION = RunParameters[i].getRunDescription();
    UInt_t NEVENTS = RunParameters[i].getEvents();
    UInt_t START_EVENT = RunParameters[i].getStartEvent();
		bool DO_PEDESTALANALYSIS = RunParameters[i].doPedestalAnalysis();
		bool DO_CLUSTERANALYSIS  = RunParameters[i].doClusterAnalysis();
		bool DO_SELECTIONANALYSIS = RunParameters[i].doSelectionAnalysis();
		bool DO_ALIGNMENT = RunParameters[i].doAlignment();
		bool DO_ALIGNMENTANALYSIS = RunParameters[i].doAlignmentAnalysis();
		cout << endl << endl << endl << endl;
		cout << "====================================" << endl;
		cout << "==> Starting analysis.." << endl;
		cout << "====================================" << endl << endl;
		cout << "RUNNUMBER: " << RUNNUMBER << endl;
		cout << "NEVENTS: " << NEVENTS << endl;
		cout << "RUNDESCRIPTION: " << RUNDESCRIPTION << endl;
		cout << "VERBOSITY: " << VERBOSITY << endl;
		cout << "INITIAL_EVENT: " << START_EVENT << endl;
		cout << "DO_PEDESTALANALYSIS: "<<DO_PEDESTALANALYSIS<<endl;
		cout << "DO_CLUSTERANALYSIS: "<<DO_CLUSTERANALYSIS<<endl;
		cout << "DO_SELECTIONANALYSIS: "<<DO_SELECTIONANALYSIS<<endl;
		cout << "DO_ALIGNMENT: " << DO_ALIGNMENT << endl;
    cout << "DO_ALIGNMENTANALYSIS: "<<DO_ALIGNMENTANALYSIS<<endl;
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

    TResults *currentResults =new TResults(settings);
    currentResults->Print();

		TRawEventSaver *eventSaver;
		eventSaver = new TRawEventSaver(RUNNUMBER);
		eventSaver->setSettings(settings);
		eventSaver->saveEvents(NEVENTS);
		delete eventSaver;

		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;

		sys->cd(currentDir.c_str());

			//Calculate Pedestal

			TPedestalCalculation* pedestalCalculation;
//			pedestalCalculation = new TPedestalCalculation(RUNNUMBER,NEVENTS);
			pedestalCalculation = new TPedestalCalculation(settings);
			pedestalCalculation->calculatePedestals(NEVENTS);
			pedestalCalculation->calculateSlidingPedestals(NEVENTS);
			delete pedestalCalculation;

		if(DO_PEDESTALANALYSIS){
			TAnalysisOfPedestal *analysisOfPedestal;
			analysisOfPedestal = new TAnalysisOfPedestal(settings);
			analysisOfPedestal->setResults(currentResults);
			analysisOfPedestal->doAnalysis(NEVENTS);
			delete analysisOfPedestal;
		}

		THTMLGenerator *htmlGen = new THTMLGenerator(settings);
		stringstream path;
		path<<currentDir<<"/"<<settings->getRunNumber()<<"/";
		htmlGen->setMainPath("./");//(string)(currentDir+"/16202/"));
		htmlGen->setSubdirPath("");

		htmlGen->setFileName("overview.html");
		htmlGen->setFileGeneratingPath(path.str());
		htmlGen->addSection("Pedestal","<a href=\"pedestalAnalysis/pedestal.html\">PEDESTAL</a>");
		htmlGen->addSection("Clustering","<a href=\"clustering/clustering.html\">CLUSTERING</a>");
		htmlGen->addSection("Selection","<a href=\"selections/selection.html\">SELECTION</a>");
		htmlGen->addSection("Landau","<a href=\"selectionAnalysis/landaus.html\">LANDAU-DISTRIBUTIONS</a>");
		htmlGen->generateHTMLFile();
		delete htmlGen;


//		process_mem_usage(vm2, rss2);
//		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;

		sys->cd(currentDir.c_str());
		TClustering* clustering;
		clustering=new TClustering(settings);//int seedDetSigma=10,int hitDetSigma=7,int seedDiaSigma=5, int hitDiaSigma=3);
		std::cout<<"cluster"<<endl;
		clustering->ClusterEvents(NEVENTS);
		delete clustering;



//		process_mem_usage(vm2, rss2);
//		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;
		sys->cd(currentDir.c_str());
		TSelectionClass* selectionClass;
		selectionClass=new TSelectionClass(settings);
		selectionClass->SetResults(currentResults);
		selectionClass->MakeSelection(NEVENTS);
		delete selectionClass;

		if (DO_CLUSTERANALYSIS){
			sys->cd(currentDir.c_str());
			TAnalysisOfClustering* analysisClustering;
			analysisClustering= new TAnalysisOfClustering(settings);
			analysisClustering->doAnalysis(NEVENTS);
			delete analysisClustering;
		}
		if(DO_SELECTIONANALYSIS){
		  sys->cd(currentDir.c_str());
			TAnalysisOfSelection *analysisSelection=new TAnalysisOfSelection(settings);
			analysisSelection->doAnalysis(NEVENTS);
			delete analysisSelection;
		}

		if (DO_ALIGNMENT){
			sys->cd(currentDir.c_str());
			TAlignment *alignment;
			alignment= new TAlignment(settings);
			alignment->setSettings(settings);
			//alignment->PrintEvents(1511,1501);
			process_mem_usage(vm2, rss2);
			cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;

			//alignment->createEventVectors(1000);
			process_mem_usage(vm2, rss2);
			cout << "\nMemory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;
//			alignment->setVerbosity(2);
			alignment->Align(NEVENTS);
			delete alignment;
		}

		if(DO_ALIGNMENTANALYSIS){
      sys->cd(currentDir.c_str());
      TAnalysisOfAlignment *anaAlignment;
      anaAlignment=new TAnalysisOfAlignment(settings);
      anaAlignment->doAnalysis(NEVENTS);
      delete anaAlignment;
		}
//		TTransparentAnalysis *transpAna;
//		transpAna = new TTransparentAnalysis(RUNNUMBER, *settings);
//		transpAna->analyze(NEVENTS,INITIAL_EVENT);

		currentResults->Print();
    currentResults->saveResults();
    delete currentResults;


		process_mem_usage(vm2, rss2);
		cout << "Memory usage: VM: " << vm2 << "; RSS: " << rss2 << endl;


    runWatch.Stop();
    cout<<"needed Time for Run "<<RUNNUMBER<<":"<<endl;
    runWatch.Print();
		if (settings!=NULL){
		  cout<<"delete Settings..."<<endl;
		  delete settings;
		  cout<<"DONE_SETTINGS"<<endl;
		}
	}
	cout<<"DONE with Analysis of all Runs "<<RunParameters.size()<<"from RunList.ini"<<endl;
	cout<<"time for all analysis:"<<endl;
	comulativeWatch.Print();
	cout<<"DONE_ALL"<<endl;

	return 0;
}


int ReadRunList() {
	RunInfo run;
	RunParameters.clear();
	cout << endl << "reading runlist.." << endl;
	ifstream file(runListPath.c_str());//"RunList.ini");
	if (!file) {
		cout << "An error has encountered while trying to open RunList.ini" << endl;
		return 0;
	}
	else cout << "RunList.ini" << " successfully opened." << endl << endl;

  int RunNumber;
  int Verbosity;
  int NEvents;
  UInt_t nStartEvent;
  char RunDescription[200];
  int bPedestalAnalysis;
  int bClusterAnalysis;
  int bSelectionAnalysis;
  int bAlignment;
  int bAlignmentAnalysis;
//  cout<<"start file loop"<<flush;
	while (!file.eof()) {
//	  RunDescription = "";
	  NEvents = 10000;
	  nStartEvent = 1000;
	  Verbosity=0;
		
		//get next line
		string line;
//		cout<<"getLine"<<endl;
		getline(file,line);
		
		//check if comment or empty line
		if ((line.substr(0, 1) == ";") || (line.substr(0, 1) == "#") || (line.substr(0, 1) == "/") || line.empty()) {
//		  cout<<"continue"<<endl;
			continue;
		}
//		cout<<"Read Line"<<endl;
		sscanf(line.c_str(), "%d %s %d %d %d %d %d %d %d %d", &RunNumber, RunDescription, &Verbosity, &NEvents, &nStartEvent, &bPedestalAnalysis, &bClusterAnalysis, &bSelectionAnalysis,&bAlignment,&bAlignmentAnalysis);
//		cout << "RunDescription Char: " << RunDescription[0] << endl;
		cout<<RunNumber<<endl;
		cout<<NEvents<<endl;
		cout<<nStartEvent<<":"<<bPedestalAnalysis<<bClusterAnalysis<<bSelectionAnalysis<<bAlignment<<bAlignmentAnalysis<<endl;
		run.setParameters(RunNumber,(string)RunDescription,Verbosity,NEvents,nStartEvent,bPedestalAnalysis,bClusterAnalysis,bSelectionAnalysis,bAlignment,bAlignmentAnalysis);
//		cout<<"Got new Parameters: "<<RunNumber<<endl;
		RunParameters.push_back(run);
	}
	return 1;
}
