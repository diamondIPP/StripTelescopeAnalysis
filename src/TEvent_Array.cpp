/*
 * TEvent_Array.cpp
 *
 *  Created on: 30.07.2011
 *      Author: Felix Bachmair
 */

#include "TEvent_Array.hh"

using namespace std;

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

//Aysha added class TEvent_Array_F method definitions
TEvent_Array_F::TEvent_Array_F(Int_t size)
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
TEvent_Array_F::~TEvent_Array_F() {}
Float_t TEvent_Array_F::GetChannelValue(Int_t index)
{
	return Channel_ADC_values[index];
}
Int_t TEvent_Array_F::GetSize()
{
	return Iteration_Size;
}
void TEvent_Array_F::SetChannelValue(Float_t value, Int_t index)
{
	Channel_ADC_values[index] = value;
}
void TEvent_Array_F::SetChannelValueSquared(Float_t value, Int_t index)
{
	//Channel_ADC_values_squared[index] = value*value; //Taylor, change #3 6/30, -2.4 s
	if(value==index){} //to shut
}
Float_t TEvent_Array_F::CalculateSigma(Float_t &new_mean)
{
	Float_t Size = Iteration_Size;
	//Float_t Sum = 0;
	//Float_t Sum2 = 0;
	float Sum = 0; //Taylor: change #4 7/1
	double Sum2 = 0; //Taylor
	Float_t Channel_ADC_value_i; //Taylor
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
Float_t TEvent_Array_F::CalculateSigma()
{
	Float_t Size = Iteration_Size;
	//Float_t Sum = 0;
	//Float_t Sum2 = 0;
	Float_t Sum = 0; //Taylor: change #4 7/1
	Double_t Sum2 = 0; //Taylor
	Float_t Channel_ADC_value_i; //Taylor
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
