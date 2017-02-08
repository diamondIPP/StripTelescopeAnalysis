/*
 * TRawEventReader.cpp
 *
 *  Created on: 01.11.2011
 *      Author: bachmair
 */

#include "../include/TRawEventReader.hh"


TRawEventReader::TRawEventReader(TSettings *settings){//Int_t runNumber) {
	// TODO Auto-generated constructor stub
	verbosity=0;
	EventsPerFile = 10000;
	this->run_number=settings->getRunNumber();
	current_rz_filename="";
	this->settings=settings;
}

TRawEventReader::~TRawEventReader() {
	// TODO Auto-generated destructor stub
	current_rz_file.close();
}

int TRawEventReader::ReadRawEvent(int EventNumber, bool verbose)
{
	//Declarations


	//Check that rzEvent size is consistant -- rzEvent is declared in Diamondstruct.h


	//Filename to lookup event
	std::ostringstream filename;
	//TSystem *sys =gSystem;
	//cout<<settings->getInputDir()<<endl;
	//cout<<settings->getAbsoluteInputPath()<<endl;
	filename << settings->getAbsoluteInputPath()<<"/RUN_" << run_number << "_" << EventNumber/EventsPerFile << ".rz";

	//	   cout<<filename.str()<<endl;


	//Open the desired rz file if not open already
	if(current_rz_filename!=filename.str()) {
		current_rz_file.close(); // desired filename is different so close old file
		current_rz_filename=filename.str();
		if(verbose) {
			std::cout << "**************************************" << std::endl;
			std::cout << "Opening " << filename.str() << " for read" << std::endl;
			std::cout << "**************************************" << std::endl;
		}
		current_rz_file.open(filename.str().c_str(),std::ios::in | std::ios::binary);  //The .c_str() must be added for ifstream to be able to read in the file name string.
		if (!current_rz_file) {
			std::cout<<settings->getAbsoluteInputPath()<<std::endl;
			std::cout << "File open error: " << filename.str() << " not found" << std::endl;
			return -1; //returning -1 signals Slide() to abort the pedestal calculation
		}
	}

	//Read in event data (Header, Data, and Trailer) using Event structure as specified in Diamondstuct.h
	current_rz_file.seekg(EventNumber%EventsPerFile * sizeof(rzEvent),ios::beg);
	current_rz_file.read(reinterpret_cast<char*>(&rzEvent),sizeof(rzEvent));

	//Changing Endianness of Event Header Read in Data
	uendian_swap(rzEvent.EvTrig);
	uendian_swap(rzEvent.EvNo);
	uendian_swap(rzEvent.EvPos);
	endian_swap(rzEvent.EvTag);
	endian_swap(rzEvent.EvDate);
	endian_swap(rzEvent.EvTime);
	uendian_swap(rzEvent.TrigCnt);
	uendian_swap(rzEvent.EvVmeTime);
	for (int i=0; i<8; i++)
	{
		endian_swap(rzEvent.VFasCnt[i]);
		endian_swap(rzEvent.VFasReg[i]);
	}
	endian_swap(rzEvent.EvNetTime);
	short_endian_swap(rzEvent.MeasNo);
	short_endian_swap(rzEvent.EvInMeasNo);
	endian_swap(rzEvent.Reserved[0]);
	endian_swap(rzEvent.Reserved[1]);

	//Endian Swap for Diamond Data and Outputing the test values for the data following the Telescope Reference Detectors
	for (int i=0; i<DIAMOND_MEM; i++)
		ushort_endian_swap(rzEvent.RD42[i]);

	//Swaping Endianness and then Outputing the Event Trailer data
	uendian_swap(rzEvent.Eor);

	if(false){
		cout << "\nEVENT:   "<<EventNumber<<"\nHeader dump:" << endl;
		cout << "EvTrig:    " << rzEvent.EvTrig << endl;
		cout << "EvNo:      " << rzEvent.EvNo << endl;
		cout << "EvPos:     " << rzEvent.EvPos << endl;
		cout << "EvTag:     " << rzEvent.EvTag << endl;
		cout << "EvDate:    " << std::bitset<32>(rzEvent.EvTime) << endl;
		cout << "TrigCnt:   " << rzEvent.TrigCnt << endl;
		cout << "EvVmeTime: " << std::bitset<32>(rzEvent.EvVmeTime) << endl;
		cout << "EvNetTime: " << std::bitset<32>(rzEvent.EvNetTime) << endl;
		cout << "MeasNo:    " << rzEvent.MeasNo << endl;
		cout << "EvInMeasNo:" << rzEvent.EvInMeasNo << endl;
	}
	//Reading out Event Header Data to Screen
	if(verbose) {
		cout << "Header dump:" << endl;
		cout << "EvTrig: " << rzEvent.EvTrig << endl;
		cout << "EvNo: " << rzEvent.EvNo << endl;
		cout << "EvPos: " << rzEvent.EvPos << endl;
		cout << "EvTag: " << rzEvent.EvTag << endl;
		cout << "EvDate: " << rzEvent.EvTime << endl;
		cout << "TrigCnt: " << rzEvent.TrigCnt << endl;
		cout << "EvVmeTime: " << rzEvent.EvVmeTime << endl;
		for (int j=0; j<8; j++) cout << "VFasCnt[" << j << "]: " << rzEvent.VFasCnt[j] << endl;
		for (int j=0; j<8; j++) cout << "VFasReg[" << j << "]: " << rzEvent.VFasReg[j] << endl;
		cout << "EvNetTime: " << rzEvent.EvNetTime << endl;
		cout << "MeasNo: " << rzEvent.MeasNo << endl;
		cout << "EvInMeasNo: " << rzEvent.EvInMeasNo << endl;
		cout << "Reserved[EVENT_HEADER_RESERVED_ESZ-2]: " << rzEvent.Reserved[0] << " and " << rzEvent.Reserved[1] << endl;
		cout << "Eor: " << rzEvent.Eor<< endl;
	}


	//Sorting the values of intoutput into different columns for the x and y strips of the 4 detectors

	/* ---------------------------------------------------------------------------------------------------------------------
	   As described in the original DAQ readme, the 2048 bytes for each event of the silicon telescope data is read in the
	   following order: To start, the first channel of the X layer ADC values are read in for each of the 4 detectors and
	   then so on for each channel until 256. Then the Y layer ADC values are read in a similar manner. So explicitly, channel 0
	   of the D0X is read in first, then channel 0 of D1X is next, then channel 0 of D2X, then channel 0 of D3X, then channel 1
	   of D0X and so on. This explains the the way the values are then sorted into the 256 byte arrays for each Detector layer.
	   ------------------------------------------------------------------------------------------------------------------------ */

	//store the raw adc values in class member storage
	for (int i=0; i<256; i++)
	{
		//NOTE: Realized that due to scrambling of the data (example 0x3615abcd written in header as 0x1536cdab)
		//and realizing that both silicon and diamond data is written down as 4-byte words, detectors should be

		//New way of mapping detectors
		D0X.ADC_values[i]=rzEvent.Input[4*i+3]; //The 1X and 3X detectors are actually physically flipped so we need to reverse the values with the program
		D1X.ADC_values[255-i]=rzEvent.Input[4*i+2];
		D2X.ADC_values[i]=rzEvent.Input[4*i+1];
		D3X.ADC_values[255-i]=rzEvent.Input[4*i];

		D0Y.ADC_values[i]=rzEvent.Input[4*i+3+1024];
		D1Y.ADC_values[i]=rzEvent.Input[4*i+2+1024];
		D2Y.ADC_values[i]=rzEvent.Input[4*i+1+1024];
		D3Y.ADC_values[i]=rzEvent.Input[4*i+1024];
	}
	for (int i=0; i<128;i++){
		Dia0.ADC_values[i]=rzEvent.RD42[i*2+1];
		Dia1.ADC_values[i]=rzEvent.RD42[i*2];
		//if(i<20)cout<<rzEvent.RD42[i*2]<<" "<<Dia0.ADC_values[i]<<"  ";
	}
	// cout<<endl;

	//Memory Consistancy check. If The amounof Diamond memory is not 263 bytes worth, then this will output a non sensical number.
	if (verbosity&&(EventNumber%1000 == 0))
	{
		cout << "For requested event " << EventNumber << ", the current value of EvTrig and EvPos is: " << rzEvent.EvTrig << " and " << rzEvent.EvPos << endl;
	}

	return 0; //returning 0 signals Slide() to continue with the calculation

}

TDetector_Data TRawEventReader::getD0X() const
{
	return D0X;
}

void TRawEventReader::setD0X(TDetector_Data d0X)
{
	D0X = d0X;
}

TDetector_Data TRawEventReader::getD0Y() const
{
	return D0Y;
}

TDetector_Data TRawEventReader::getD1X() const
{
	return D1X;
}

TDetector_Data TRawEventReader::getD1Y() const
{
	return D1Y;
}

TDetector_Data TRawEventReader::getD2X() const
{
	return D2X;
}

TDetector_Data TRawEventReader::getD2Y() const
{
	return D2Y;
}

TDetector_Data TRawEventReader::getD3X() const
{
	return D3X;
}

TDetector_Data TRawEventReader::getD3Y() const
{
	return D3Y;
}

TDetector_Data TRawEventReader::getDia0() const
{
	return Dia0;
}

TDetector_Data TRawEventReader::getDia1() const
{
	return Dia1;
}

void TRawEventReader::setD0Y(TDetector_Data d0Y)
{
	D0Y = d0Y;
}

void TRawEventReader::setD1X(TDetector_Data d1X)
{
	D1X = d1X;
}

void TRawEventReader::setD1Y(TDetector_Data d1Y)
{
	D1Y = d1Y;
}

void TRawEventReader::setD2X(TDetector_Data d2X)
{
	D2X = d2X;
}

void TRawEventReader::setD2Y(TDetector_Data d2Y)
{
	D2Y = d2Y;
}

void TRawEventReader::setD3X(TDetector_Data d3X)
{
	D3X = d3X;
}

void TRawEventReader::setD3Y(TDetector_Data d3Y)
{
	D3Y = d3Y;
}

void TRawEventReader::setDia0(TDetector_Data dia0)
{
	Dia0 = dia0;
}

TDetector_Data TRawEventReader::getPlane(int det,UInt_t diaInput)
{

	switch(det){
	case 0: return D0X;break;
	case 1: return D0Y;break;
	case 2: return D1X;break;
	case 3: return D1Y;break;
	case 4: return D2X;break;
	case 5: return D2Y;break;
	case 6: return D3X;break;
	case 7: return D3Y;break;
	case 8: if(diaInput==1)
		return Dia1;
	else
		return Dia0;
	break;
	}
	cerr<<" Trying to get data from detector "<<det<<", which does not exist!!!!!"<<endl;
	return TDetector_Data();
}

TDetector_Data TRawEventReader::getDia(UInt_t diaInput) const
{
	if(diaInput==0)
		return getDia0();
	else
		return getDia1();
}

void TRawEventReader::setDia1(TDetector_Data dia1)
{
	Dia1 = dia1;
}



//end of binary data read-in while loop


