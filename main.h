int RUNNUMBER;
string RUNDESCRIPTION;
int NEVENTS;
int INITIAL_EVENT;
int HIT_OCCUPANCY;
bool PLOTS;
bool ALTERNATIVECLUSTERING;
bool RunListOK;

void initVariables();
int ReadRunList();

class RunInfo {
public:
	int RunNumber, NEvents, Initial_Event, Hit_Occupancy;
	string RunDescription;
	bool AlternativeClustering;
	void SetParameters() {
		RunNumber = RUNNUMBER;
		NEvents = NEVENTS;
		Initial_Event = INITIAL_EVENT;
		Hit_Occupancy = HIT_OCCUPANCY;
		RunDescription = RUNDESCRIPTION;
		AlternativeClustering = ALTERNATIVECLUSTERING;
	}
	void GetParameters() {
		RUNNUMBER = RunNumber;
		NEVENTS = NEvents;
		INITIAL_EVENT = Initial_Event;
		HIT_OCCUPANCY = Hit_Occupancy;
		RUNDESCRIPTION = RunDescription;
		ALTERNATIVECLUSTERING = AlternativeClustering;
	}
};
vector<RunInfo> RunParameters;