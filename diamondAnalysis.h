#include <fstream>
#include <iostream>
int RUNNUMBER;
std::string RUNDESCRIPTION;
int NEVENTS;
int INITIAL_EVENT;
int HIT_OCCUPANCY;
int VERBOSITY;
bool PLOTS;
bool ALTCLUSTERING;
bool RunListOK;
int DO_SLIDINGPEDESTAL, DO_ALIGNMENT, CUTFAKETRACKS;
std::string inputDir;
std::string outputDir;
std::string runListPath="RunList.ini";

void initVariables();
int ReadRunList();

class RunInfo {
public:
	int RunNumber, Verbosity, NEvents, Initial_Event, Hit_Occupancy;
	std::string RunDescription;
	bool AlternativeClustering, DoAlignment, DoSlidingPedestal, CutFakeTracks;
	void SetParameters() {
		RunNumber = RUNNUMBER;
        Verbosity = VERBOSITY;
		NEvents = NEVENTS;
		Initial_Event = INITIAL_EVENT;
		Hit_Occupancy = HIT_OCCUPANCY;
		RunDescription = RUNDESCRIPTION;
		AlternativeClustering = ALTCLUSTERING;
		DoAlignment = DO_ALIGNMENT;
		DoSlidingPedestal = DO_SLIDINGPEDESTAL;
		CutFakeTracks = CUTFAKETRACKS;
	}
	void GetParameters() {
		RUNNUMBER = RunNumber;
        VERBOSITY = Verbosity;
		NEVENTS = NEvents;
		INITIAL_EVENT = Initial_Event;
		HIT_OCCUPANCY = Hit_Occupancy;
		RUNDESCRIPTION = RunDescription;
		ALTCLUSTERING = AlternativeClustering;
		DO_ALIGNMENT = DoAlignment;
		DO_SLIDINGPEDESTAL = DoSlidingPedestal;
		CUTFAKETRACKS = CutFakeTracks;
	}
};
std::vector<RunInfo> RunParameters;
