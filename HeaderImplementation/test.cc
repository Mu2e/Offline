#include "FindMultiplePeak.hh"
#include <iostream>
#include "FitModelRoot.hh" 
#include "TF1.h"
#include "FindPeakBase.hh"

void testit()
{
	adcWaveform adcData;
	unsigned int adc[8] = {64, 60, 103, 125, 116, 108, 93, 92};
	adcData = adc;
	resultantHitData result;
	const configStruct init;
	FindMultiplePeaks peak(init);
	peak.process(adcData, result);
	std::cout << result[0]._peakTime << std::endl;
	std::cout << result[0]._peakHeight << std::endl;
}

void testit2(Double_t timeShift)
{
	TF1 *func = new TF1("func",FitModelRoot::convolutionSinglePeakWithConstantPedestal,0.0,140.0,4);
	func->SetParameters(timeShift,8494,0.0,10.0);
	Double_t adc[8] = {64, 60, 103, 125, 116, 108, 93, 92};
	Double_t times[8] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0};
	TGraph * gr = new TGraph(8, times, adc);
	gr->Fit(func);
	gr->Draw("A*");
	func->Draw("same");
}


