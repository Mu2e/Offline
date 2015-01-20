{
	#include "classImplementation.C"
	convolutionSinglePeakWithDynamicPedestal *c = new convolutionSinglePeakWithDynamicPedestal();
	unsigned int adc[8] = {64, 60, 103, 125, 116, 108, 93, 92};

	paramStruct initParams;


	FindSinglePeak *d = new FindSinglePeak(initParams, adc);
	d->adcPeaks.push_back(125.0);
	d->timePeaks.push_back(60.0);

	//d->process();

}