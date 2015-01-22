#include "classImplementation.C"
#include <iostream>

void testit()
{
	unsigned int adc[8] = {64, 60, 103, 125, 116, 108, 93, 92};

	configStruct initParams;

	resultantHitData result;
//	resultantPeakData peakData(60.0,125.0);
//	result.push_back(peakData);
	FindMultiplePeaks d(initParams);
	d.process(result,adc);
	std::cout << result[0].peakTime << std::endl;
	std::cout << result[0].peakHeight << std::endl;
}