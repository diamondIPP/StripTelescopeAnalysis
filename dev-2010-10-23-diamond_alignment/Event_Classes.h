//2010-07-12 Taylor's tweaks added
//2010-08-15 Taylor fixed problem needed long instead of int
//2010-09-18 Changed Sum2 from int to long

#ifndef EventClasses
#define EventClasses

#include "TMath.h"
#include <iostream>

using namespace std;

class TDetector_Data {

   public: 
      TDetector_Data();
      ~TDetector_Data();
      int GetADC_value(Int_t index);
      void SetADC_value(Int_t index, Int_t value);
      int ADC_values[256];

};

TDetector_Data::TDetector_Data()
{
   for(Int_t i=0; i<256; i++) 
   {
      ADC_values[i]=0;
   }
}

TDetector_Data::~TDetector_Data() {}

int TDetector_Data::GetADC_value(Int_t index)
{
   return ADC_values[index];
}

void TDetector_Data::SetADC_value(Int_t index, Int_t value)
{
   ADC_values[index] = value;
}

class TTrigger_Event {

   public:
      TTrigger_Event();
      ~TTrigger_Event();
      //TDetector_Data GetDetector(TDetector_Data Detector);
      void SetD0X(TDetector_Data Detector);
      void SetD0Y(TDetector_Data Detector);
      void SetD1X(TDetector_Data Detector);
      void SetD1Y(TDetector_Data Detector);
      void SetD2X(TDetector_Data Detector);
      void SetD2Y(TDetector_Data Detector);
      void SetD3X(TDetector_Data Detector);
      void SetD3Y(TDetector_Data Detector);
      void SetDia0(TDetector_Data Detector);
      void SetDia1(TDetector_Data Detector);
      void SetAny(TDetector_Data Detector);
      TDetector_Data GetD0X() const {return D0X;}
      TDetector_Data GetD0Y() const {return D0Y;}
      TDetector_Data GetD1X() const {return D1X;}
      TDetector_Data GetD1Y() const {return D1Y;}
      TDetector_Data GetD2X() const {return D2X;}
      TDetector_Data GetD2Y() const {return D2Y;}
      TDetector_Data GetD3X() const {return D3X;}
      TDetector_Data GetD3Y() const {return D3Y;}
      TDetector_Data GetDia0() const {return Dia0;}
      TDetector_Data GetDia1() const {return Dia1;}
      TDetector_Data GetAny() const {return Any;}

   protected:
      TDetector_Data D0X;
      TDetector_Data D0Y;
      TDetector_Data D1X;
      TDetector_Data D1Y;
      TDetector_Data D2X;
      TDetector_Data D2Y;
      TDetector_Data D3X;
      TDetector_Data D3Y;
      TDetector_Data Dia0;
      TDetector_Data Dia1;
      TDetector_Data Any;

};

TTrigger_Event::TTrigger_Event() {}
TTrigger_Event::~TTrigger_Event() {}
void TTrigger_Event::SetD0X(TDetector_Data Detector)
{
   D0X = Detector;
}
void TTrigger_Event::SetD0Y(TDetector_Data Detector)
{
   D0Y = Detector;
}
void TTrigger_Event::SetD1X(TDetector_Data Detector)
{
   D1X = Detector;
}
void TTrigger_Event::SetD1Y(TDetector_Data Detector)
{
   D1Y = Detector;
}
void TTrigger_Event::SetD2X(TDetector_Data Detector)
{
   D2X = Detector;
}
void TTrigger_Event::SetD2Y(TDetector_Data Detector)
{
   D2Y = Detector;
}
void TTrigger_Event::SetD3X(TDetector_Data Detector)
{
   D3X = Detector;
}
void TTrigger_Event::SetD3Y(TDetector_Data Detector)
{
   D3Y = Detector;
}
void TTrigger_Event::SetDia0(TDetector_Data Detector)
{
   Dia0 = Detector;
}
void TTrigger_Event::SetDia1(TDetector_Data Detector)
{
   Dia1 = Detector;
}
void TTrigger_Event::SetAny(TDetector_Data Detector)
{
   Any = Detector;
}




class TPed_and_RMS {

   public:
      TPed_and_RMS();
      ~TPed_and_RMS();
      void SetPedValues(Int_t ped_index, Float_t ped_value);
      void SetRMSValues(Int_t RMS_index, Float_t RMS_value);
      Float_t GetPedValues(Int_t ped_index){return Pedestal_values[ped_index];}
      Float_t GetRMSValues(Int_t RMS_index){return RMS_values[RMS_index];}

   protected:
      Float_t Pedestal_values[256];
      Float_t RMS_values[256];
};
TPed_and_RMS::TPed_and_RMS()
{
   for(Int_t i=0; i<256; i++)
   {
      Pedestal_values[i] = 0;
      RMS_values[i] = 1000;
   }
}
TPed_and_RMS::~TPed_and_RMS() {}
void TPed_and_RMS::SetPedValues(Int_t ped_index, Float_t ped_value)
{
   Pedestal_values[ped_index]=ped_value;
}
void TPed_and_RMS::SetRMSValues(Int_t RMS_index, Float_t RMS_value)
{
   RMS_values[RMS_index]=RMS_value;
}

class TEvent_Array {
   public:
      TEvent_Array(Int_t size);
      ~TEvent_Array();
      Int_t GetChannelValue(Int_t index);
      Int_t GetSize();
      void SetChannelValue(Int_t value, Int_t index);
      void SetChannelValueSquared(Int_t value, Int_t index);
      Float_t CalculateSigma(Float_t &new_mean);
      Float_t CalculateSigma();

   private:
      Int_t Iteration_Size;
      vector<Int_t> Channel_ADC_values;
      //Int_t Channel_ADC_values[1000];
      //Int_t Channel_ADC_values_squared[1000]; //Taylor
};
TEvent_Array::TEvent_Array(Int_t size)
{
  Iteration_Size = size;
  //for(Int_t i=0; i<300; i++)
  for (Int_t i = 0; i < Iteration_Size; i++)
  {
      Channel_ADC_values.push_back(0);
      //Channel_ADC_values[i]=0;
      //Channel_ADC_values_squared[i]=0; //Taylor
   }
}
TEvent_Array::~TEvent_Array() {}
Int_t TEvent_Array::GetChannelValue(Int_t index)
{
   return Channel_ADC_values[index];
}
Int_t TEvent_Array::GetSize()
{
   return Iteration_Size;
}
void TEvent_Array::SetChannelValue(Int_t value, Int_t index)
{
   Channel_ADC_values[index] = value;
}
void TEvent_Array::SetChannelValueSquared(Int_t value, Int_t index)
{
   //Channel_ADC_values_squared[index] = value*value; //Taylor, change #3 6/30, -2.4 s
   if(value==index){} //to shut 
}
Float_t TEvent_Array::CalculateSigma(Float_t &new_mean)
{
   Float_t Size = Iteration_Size;
   //Float_t Sum = 0;
   //Float_t Sum2 = 0;
   int Sum = 0; //Taylor: change #4 7/1
   long Sum2 = 0; //Taylor
   Int_t Channel_ADC_value_i; //Taylor
   for(Int_t i=0; i<Size; i++)
   {
      //Sum = Sum+Channel_ADC_values[i];
      //Sum2 = Sum2+Channel_ADC_values_squared[i];
      Channel_ADC_value_i = Channel_ADC_values[i];
      Sum += Channel_ADC_value_i;
      Sum2 += Channel_ADC_value_i * Channel_ADC_value_i;

      if (Sum < 0) cout << "TEvent_Array: Warning!! Mistake: Sum = " << Sum << endl; //Taylor
      if (Sum2 < 0) cout << "TEvent_Array: Warning!! Mistake: Sum2 = " << Sum2 << endl; //Taylor
   }
   new_mean = Sum / Size;
   //return TMath::Sqrt((Sum2/Size)-((Sum/Size)*(Sum/Size)));
   return TMath::Sqrt(Sum2 / Size - new_mean * new_mean);
   //return TMath::Sqrt((Sum2 - Sum * Sum / Size) / Size); //Taylor - gave a different result -> indicated lack of precision

   /*Float_t option1 = (Sum2/Size)-((Sum/Size)*(Sum/Size));
   Float_t option2 = (Sum2 - Sum * Sum / Size) / Size;
   if (TMath::Abs(option1 - option2) > 1)
   cout << option1 << " : " << option2 << endl; // */ //Taylor - outputs nothing

   //Taylor: this somehow ends up making a big difference in the final noise calculation
   //Solution: Sum/2 float -> int
}
Float_t TEvent_Array::CalculateSigma()
{
   Float_t Size = Iteration_Size;
   //Float_t Sum = 0;
   //Float_t Sum2 = 0;
   Int_t Sum = 0; //Taylor: change #4 7/1
   Long64_t Sum2 = 0; //Taylor
   Int_t Channel_ADC_value_i; //Taylor
   for(Int_t i=0; i<Size; i++)
   {
      //Sum = Sum+Channel_ADC_values[i];
      //Sum2 = Sum2+Channel_ADC_values_squared[i];
      Channel_ADC_value_i = Channel_ADC_values[i];
      Sum += Channel_ADC_value_i;
      Sum2 += Channel_ADC_value_i * Channel_ADC_value_i;

      if (Sum < 0) cout << "TEvent_Array: Warning!! Mistake: Sum = " << Sum << endl; //Taylor
      if (Sum2 < 0) cout << "TEvent_Array: Warning!! Mistake: Sum2 = " << Sum2 << endl; //Taylor
   }
   return TMath::Sqrt((Sum2/Size)-((Sum/Size)*(Sum/Size)));
   //return TMath::Sqrt((Sum2 - Sum * Sum / Size) / Size); //Taylor - gave a different result -> indicated lack of precision

   /*Float_t option1 = (Sum2/Size)-((Sum/Size)*(Sum/Size));
   Float_t option2 = (Sum2 - Sum * Sum / Size) / Size;
   if (TMath::Abs(option1 - option2) > 1)
   cout << option1 << " : " << option2 << endl; // */ //Taylor - outputs nothing

   //Taylor: this somehow ends up making a big difference in the final noise calculation
   //Solution: Sum/2 float -> int
}


#endif
