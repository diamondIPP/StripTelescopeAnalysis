int RUNNUMBER;
string RUNDESCRIPTION;
int NEVENTS;
int INITIAL_EVENT;
int HIT_OCCUPANCY;
bool PLOTS;
bool ALTERNATIVECLUSTERING;
bool RunListOK;
bool DO_SLIDINGPEDESTAL, DO_ALIGNMENT;

void initVariables();
int ReadRunList();

class RunInfo {
public:
	int RunNumber, NEvents, Initial_Event, Hit_Occupancy;
	string RunDescription;
	bool AlternativeClustering, DoAlignment, DoSlidingPedestal;
	void SetParameters() {
		RunNumber = RUNNUMBER;
		NEvents = NEVENTS;
		Initial_Event = INITIAL_EVENT;
		Hit_Occupancy = HIT_OCCUPANCY;
		RunDescription = RUNDESCRIPTION;
		AlternativeClustering = ALTERNATIVECLUSTERING;
		DoAlignment = DO_ALIGNMENT;
		DoSlidingPedestal = DO_SLIDINGPEDESTAL;
	}
	void GetParameters() {
		RUNNUMBER = RunNumber;
		NEVENTS = NEvents;
		INITIAL_EVENT = Initial_Event;
		HIT_OCCUPANCY = Hit_Occupancy;
		RUNDESCRIPTION = RunDescription;
		ALTERNATIVECLUSTERING = AlternativeClustering;
		DO_ALIGNMENT = DoAlignment;
		DO_SLIDINGPEDESTAL = DoSlidingPedestal;
	}
};
vector<RunInfo> RunParameters;