{
	#include "classImplementation.C"
	unsigned int adc[8] = {64, 60, 103, 125, 116, 108, 93, 92};

	configStruct initParams;

	resultantHitData result;
	resultantPeakData peakData(60.0,125.0);
	result.push_back(peakData);
	FindSinglePeak d(initParams);
	d.process(result,adc);

}