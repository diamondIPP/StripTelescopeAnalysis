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
#include "TADCEventReader.hh"
#include "TPedestalCalculation.hh"
#include "TAnalysisOfPedestal.hh"
#include "TAnalysisOfClustering.hh"
#include "TAnalysisOfSelection.hh"
#include "TAnalysisOf3dDiamonds.hh"
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
#include "TRunInfo.hh"
#include "TFile.h"
#include <string>     // std::string, std::stoi
#include <stdlib.h>     /* strtol */

using namespace std;
/*** USAGE ***/
void printHelp( void )
{
	cout<<"*******************************************************"<<endl;
	cout<<"diamondAnalysis: A Tool for anaylsis of RD42 test beams"<<endl;
	cout<<"USEAGE:\n"<<endl;
	cout<<"diamondAnalysis "<<endl;
	cout<<" Options:"<<endl;
	cout<<"\t-i INPUTDIR \n\t-o OUTPUTDIR\n\t-s RunSettingsDIR"<<endl;
	cout<<"*******************************************************"<<endl;
}

bool checkDir(string dir){
	struct stat st;
	if(stat(dir.c_str(),&st) == 0){
		printf(" %s is present\n",dir.c_str());
		return true;
	}
	cout<<dir << "seems not to be present. What to do?Exit? [y/n]\t"<<flush;
	char t='l';
	cin >>t;
	string possibleInputs = "yYnN";
	while (possibleInputs.find(t)==string::npos){
		cout<<"Please enter a valid input: 'y,Y,n,N'"<<flush;
		cin >> t;
	}
	if (t=='y'||t=='Y')
		exit(-1);
	else return false;

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
			cout<<"found outputDir: \""<<outputDir<<"\""<<endl;
		}
		if((string(argv[i]) == "-r"||string(argv[i])=="-R")&&i+1<argc){
			i++;
			runListPath=string(argv[i]);
			cout<<"runListpath is set to:\""<<runListPath<<"\""<<endl;
		}

		if((string(argv[i]) == "-s"||string(argv[i])== "-S")&& i+1<argc){
			i++;
			runSettingsDir=string(argv[i]);
			cout<<"settingDirPath is set to: \""<<runSettingsDir<<"\""<<endl;
		}
        if ( ( (string(argv[i]) == "-3d") || (string(argv[i]) == "-3D")) && i+1 < argc){
                i++;
                //run_3danalysis=std::stoi(argv[i]);
                run_3danalysis= strtol(string(argv[i]).c_str(),0,10);
                cout<<"Running 3D Analysis: "<<run_3danalysis<<endl;
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
	if(!RunListOK)exit(-1);
	TSystem* sys = gSystem;
	std::string currentDir = sys->pwd();
	for (unsigned int i = 0; i < RunParameters.size(); i++) {
		cout << RunParameters[i].getRunNumber();
		if (i+1 < RunParameters.size()) cout << ", ";
	}
	cout << " will be analysed.." << endl;
	int bufsize = 300;
	char buf[bufsize];
	readlink("/proc/self/exe", buf, bufsize);
	TString path = buf;
	path = path(0,path.Last('/'));

	if (!RunListOK) return 0;

	/**Start with Analyising, read RunParameteres of the Run and start analysis with that parameters
	 */
	for (unsigned int i = 0; i < RunParameters.size(); i++) {
		TStopwatch runWatch;
		runWatch.Start(true);
		RunParameters[i].Print();
		bool DO_ALIGNMENT = RunParameters[i].doAlignment();
		bool DO_ALIGNMENTANALYSIS = RunParameters[i].doAlignmentAnalysis();
		bool DO_TRANSPARENT_ANALYSIS = RunParameters[i].doTransparentAnalysis();

		time_t rawtime;
		tm *timestamp;
		time (&rawtime);
		timestamp = gmtime(&rawtime);

		ostringstream logfilename;
		logfilename << "analyse_log_" << RunParameters[i].getRunNumber() << "_" << timestamp->tm_year << "-" << timestamp->tm_mon << "-" << timestamp->tm_mday << "." << timestamp->tm_hour << "." << timestamp->tm_min << "." << timestamp->tm_sec << ".log";
		//
		//		FILE *log;
		//		log = freopen(logfilename.str().c_str(), "w", stdout);


		TSettings *settings = 0;
		cout<<"settings"<<endl;
		settings = new TSettings((TRunInfo*)&RunParameters[i]);
        cout<<settings->diamondPattern.hasInvalidIntervals()<<endl;

		TResults *currentResults =new TResults(settings);
		currentResults->Print();


		TRawEventSaver *eventSaver;
		eventSaver = new TRawEventSaver(settings);
		eventSaver->saveEvents(RunParameters[i].getEvents());
		delete eventSaver;

		Long64_t skipped = 0;
		if(settings->isEventSkip()){
			for(int ii = 0; ii < settings->getSkipEvents().size(); ii++){
				skipped += (settings->getSkipEvents().at(ii).second - settings->getSkipEvents().at(ii).first);
			}
		}
		Long64_t nEvents = RunParameters[i].getEvents() - skipped; // DA: +1 ?

		//Calculate Pedestal
		sys->cd(currentDir.c_str());
		TPedestalCalculation* pedestalCalculation;
		pedestalCalculation = new TPedestalCalculation(settings);
		pedestalCalculation->calculateSlidingPedestals(nEvents);
		delete pedestalCalculation;

		if(RunParameters[i].doPedestalAnalysis()){
			sys->cd(currentDir.c_str());
			TAnalysisOfPedestal *analysisOfPedestal;
			analysisOfPedestal = new TAnalysisOfPedestal(settings);
			analysisOfPedestal->setResults(currentResults);
			analysisOfPedestal->doAnalysis(nEvents);
			delete analysisOfPedestal;
		}

		THTMLGenerator *htmlGen = new THTMLGenerator(settings);
		htmlGen->setFileGeneratingPath(settings->getAbsoluteOuputPath(true));
		htmlGen->setMainPath("./");
		htmlGen->setSubdirPath("");
		htmlGen->setFileName("overview.html");
		htmlGen->addSection("Pedestal","<a href=\"./pedestalAnalysis/pedestal.html\">PEDESTAL</a>");
		htmlGen->addSection("Clustering","<a href=\"./clustering/clustering.html\">CLUSTERING</a>");
		htmlGen->addSection("Selection","<a href=\"./selections/selection.html\">SELECTION</a>");
		htmlGen->addSection("Alignment","<a href=\"./alignment/alignment.html\">ALIGNMENT</a>");
		htmlGen->addSection("Landau","<a href=\"./selectionAnalysis/landaus.html\">LANDAU-DISTRIBUTIONS</a>");
		htmlGen->generateHTMLFile();
		delete htmlGen;

		sys->cd(currentDir.c_str());
		TClustering* clustering;
		clustering=new TClustering(settings);//int seedDetSigma=10,int hitDetSigma=7,int seedDiaSigma=5, int hitDiaSigma=3);
		clustering->ClusterEvents(nEvents);
		delete clustering;

		if(RunParameters[i].doClusterAnalysis()){
			sys->cd(currentDir.c_str());
			TAnalysisOfClustering* analysisClustering;
			analysisClustering = new TAnalysisOfClustering(settings);
			analysisClustering->setResults(currentResults);
			analysisClustering->doAnalysis(nEvents);
			delete analysisClustering;
		}

		sys->cd(currentDir.c_str());
		TSelectionClass* selectionClass = new TSelectionClass(settings);
		selectionClass->SetResults(currentResults);
		selectionClass->MakeSelection(nEvents);
		delete selectionClass;
		currentResults->createResultFiles();

		if(RunParameters[i].doSelectionAnalysis()){
			sys->cd(currentDir.c_str());
			TAnalysisOfSelection *analysisSelection=new TAnalysisOfSelection(settings);
			analysisSelection->doAnalysis(nEvents);
			delete analysisSelection;
			currentResults->createResultFiles();
		}

		if (DO_ALIGNMENT){
			sys->cd(currentDir.c_str());
			TAlignment *alignment = new TAlignment(settings);
			//			alignment->setSettings(settings);
			//alignment->PrintEvents(1511,1501);
			alignment->Align(nEvents,0,TAlignment::enumDetectorsToAlign(settings->getAlignmentMode()));
			delete alignment;
			currentResults->createResultFiles();
		}
		//		if(settings->is3dDiamond()){
        //

        if ( run_3danalysis ){
            if(settings->do3dShortAnalysis()||settings->do3dLongAnalysis()||settings->do3dTransparentAnalysis()){
                TAnalysisOf3dDiamonds* analyse3dDiamond = new TAnalysisOf3dDiamonds(settings);
                analyse3dDiamond->doAnalysis(nEvents);
                delete analyse3dDiamond;
                if (settings->do3dTransparentAnalysis())
                    settings->set3dTransparentAnalysis(0);
                else
                    settings->set3dTransparentAnalysis(1);
                analyse3dDiamond = new TAnalysisOf3dDiamonds(settings);
                analyse3dDiamond->doAnalysis(nEvents);
                delete analyse3dDiamond;
                currentResults->createResultFiles();
            }
        }

		if(DO_ALIGNMENTANALYSIS){
			sys->cd(currentDir.c_str());
			TAnalysisOfAlignment *anaAlignment;
			anaAlignment=new TAnalysisOfAlignment(settings);
			anaAlignment->doAnalysis(nEvents);
			delete anaAlignment;
			currentResults->createResultFiles();
		}

		if (DO_TRANSPARENT_ANALYSIS) {
			TTransparentAnalysis *transpAna;
			transpAna = new TTransparentAnalysis(settings);
			transpAna->setResults(currentResults);
			transpAna->analyze(nEvents,RunParameters.at(i).getStartEvent());
			delete transpAna;
			currentResults->createResultFiles();
		}

		if ((DO_TRANSPARENT_ANALYSIS||DO_ALIGNMENT) && settings && settings->doTransparentAlignment()){
			sys->cd(currentDir.c_str());
			TAlignment *alignment = new TAlignment(settings,TSettings::transparentMode);
			alignment->createTransparentEventVectors(nEvents);
			//			alignment->setSettings(settings);
			//alignment->PrintEvents(1511,1501);
			alignment->Align(nEvents,0,TAlignment::diaAlignment);
			delete alignment;
			currentResults->createResultFiles();

			TTransparentAnalysis *transpAna;
			transpAna = new TTransparentAnalysis(settings,TSettings::transparentMode);
			transpAna->setResults(currentResults);
			transpAna->analyze(nEvents,RunParameters.at(i).getStartEvent());
			delete transpAna;
			currentResults->createResultFiles();
		}
		cout<<"PRINT RESULTS"<<endl;
		currentResults->Print();
		cout<<"SAVE RESULTS"<<endl;
//		currentResults->saveResults(settings->getResultsRootFilePath());
		TFile* file = new TFile (settings->getResultsRootFilePath().c_str(),"RECREATE");
		file->Print();
		file->ls();
		file->SetName("fileName");
		file->cd();
		cout<<"KEYS: "<<file->GetNkeys()<<endl;
		file->GetListOfKeys()->Print();
		cout<<"Write"<<endl;

		currentResults->Write("test");
		cout<<"CLOSE FILE"<<endl;
		file->Close();
//		if (currentResults){
//			cout<<"DELETE RESULTS"<<endl;
//			delete currentResults;
//			cout<<"#"<<endl;;
//		}
		cout<<"saved results..."<<endl;


		runWatch.Stop();
		cout<<"needed Time for Run "<<RunParameters[i].getRunNumber()<<":"<<endl;
		runWatch.Print();
   
        cout<<"needed Time for Run "<<RunParameters[i].getRunNumber()<<":\n\t"<<flush;
        runWatch.Print();
        if (settings!=NULL){
            cout<<"delete Settings..."<<endl;
            delete settings;
            cout<<"DONE_SETTINGS"<<endl;
                                                           
        }
        cout<<"DONE with Analysis of Run "<< RunParameters[i].getRunNumber();
        cout<<": "<<i+1<<"/"<<RunParameters.size()<<endl;
        cout<<RunParameters[i].getOutputDir()<<endl;
        cout<<"Bla"<<std::endl;

        //settings->getAbsoluteOuputPath(true)<<"\""<<std::endl;
        //find /scratch/analysis/output/1* -type d -exec cp -v /scratch/analysis/output/index.php {}
        TString out_path = RunParameters[i].getOutputDir();
        cout<<out_path<<endl;
        out_path =out_path+TString::Format("/%d",RunParameters[i].nRunNumber);
        cout<<out_path<<endl;
        cout<<endl;
        TString cmd = TString("find ")+out_path+ TString(" -type d -exec cp ")+path+TString("/index.php {} \\;");
        cout<<"Copying indexfile into all subdirectories: "<<cmd<<endl;
        system(cmd);
   }
   cout<<"DONE with Analysis of all "<<RunParameters.size()<<" Runs from RunList.ini"<<endl;
   cout<<"total time for all analysis:\n\t"<<flush;

	cout<<"time for all analysis:"<<endl;
	comulativeWatch.Print();
	cout<<"DONE_ALL. Run with Revision: "<<SVN_REV<<endl;

	return 0;
}


int ReadRunList() {
	TRunInfo run;
	RunParameters.clear();
	cout << endl << "reading runlist.."  << flush;
	cout << "RunListPath: \""<<runListPath<<"\""<<endl;

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
	int bTransAna;
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
		sscanf(line.c_str(), "%d %s %d %d %d %d %d %d %d %d %d", &RunNumber, RunDescription, &Verbosity, &NEvents, &nStartEvent, &bPedestalAnalysis, &bClusterAnalysis, &bSelectionAnalysis,&bAlignment,&bAlignmentAnalysis, &bTransAna);
		//		cout << "RunDescription Char: " << RunDescription[0] << endl;
		cout<<RunNumber<<endl;
		cout<<NEvents<<endl;
		cout<<nStartEvent<<":"<<bPedestalAnalysis<<bClusterAnalysis<<bSelectionAnalysis<<bAlignment<<bAlignmentAnalysis<<endl;
		run.setParameters(RunNumber,(string)RunDescription,Verbosity,NEvents,nStartEvent,bPedestalAnalysis,bClusterAnalysis,bSelectionAnalysis,bAlignment,bAlignmentAnalysis,bTransAna);
		run.setRunSettingsDir(runSettingsDir);
		run.setOutputDir(outputDir);
		run.setInputDir(inputDir);
		cout<<"output dir: "<<run.getOutputDir()<<endl;
		cout<<"input dir: "<<run.getInputDir()<<endl;
		RunParameters.push_back(run);
	}
	return 1;
}
