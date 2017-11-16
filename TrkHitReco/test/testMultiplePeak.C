#include "Mu2e3/inc/FindMultiplePeak.hh"
#include <iostream>
#include "Mu2e3/inc/ConfigStruct.hh"

const ConfigStruct initParams;

// Testing for implementation of SinglePeak class

// Processes sample data with fixed pedestal

// Obtained from file "electronFit25ns_8UniformFirst.root"
// Hit 6
bool testSinglePeak()
{
	// Create test waveform 
	adcWaveform adcData;
//	unsigned int adc[8] = {60, 64, 69, 208, 252, 198, 152, 111};

	unsigned int adc[8] = {64, 60, 103, 125, 116, 108, 93, 92};
	adcData = adc;

	resultantHitData result;
	FindSinglePeak peak(initParams);
	std::cout << "This code is working" << std::endl;
	Double_t fitParameters[5] = {1.0,2.0,3.0,4.0,5.0};
	peak.fitParams2ResultantData(fitParameters,result);
//	peak.process(adcData, result);
	// peakTime = 48.7063
	// peakHeight = 191.584
    //std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    //std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;

    return true;
}
/**
// Hit 35
bool testSinglePeakDynamicPedestal()
{
	// Create test waveform 
	adcWaveform adcData;
	unsigned int adc[8] = {124, 92, 139, 123, 107, 93, 79, 75};
	adcData = adc;

	resultantHitData result;

	FindSinglePeakWithDynamicPedestal peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;

    return true;
}


// Hit 165
bool testSinglePeakConstantPedestal()
{
	// Create test waveform 
	adcWaveform adcData;
	unsigned int adc[8] = {120, 114, 116, 241, 274, 225, 165, 120};
	adcData =  adc;

	resultantHitData result;

	FindSinglePeakWithConstantPedestal peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;

    return true;
}



// Hit 134
bool testDoublePeak()
{
	// Create test waveform 
	adcWaveform adcData;
	unsigned int adc[8] = {67, 64, 130, 189, 159, 126, 129, 599};
	adcData = adc;

	resultantHitData result;

	FindDoublePeak peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;

    return true;
}

// Hit 2454
bool testDoublePeakDynamicPedestal()
{
	// Create test waveform 
	adcWaveform adcData;
	unsigned int adc[8] = {131, 99, 102, 234, 266, 210, 241, 365};
	adcData = adc;

	resultantHitData result;

	FindDoublePeakWithDynamicPedestal peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;
 	std::cout << "result[0]._peakTime" << result[2]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[2]._peakHeight << std::endl;

    return true;
}

// Hit 2351
// Note that this case is EXTREMELY rare (~3 / 10000 hits) with no good results
bool testDoublePeakConstantPedestal()
{
	// Create test waveform 
	adcWaveform adcData;
	unsigned int adc[8] = {80, 72, 213, 388, 350, 536, 530, 384};
	adcData = adc;

	resultantHitData result;

	FindDoublePeakWithDynamicPedestal peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;
 	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;

    return true;

}


// Testing for implementation of Multiple Peaks

// Obtained from file "electronFit25ns_8UniformFirst.root"
// Hit 6
bool searchForSinglePeak()
{
	adcWaveform adcData;
	unsigned int adc[8] = {60, 64, 69, 208, 252, 198, 152, 111};
	adcData = adc;

	resultantHitData result;
	FindMultiplePeaks peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;

    return true;
}


bool searchForSinglePeakConstantPedestal()
{
	adcWaveform adcData;
	unsigned int adc[8] = {120, 114, 116, 241, 274, 225, 165, 120};
	adcData = adc;

	resultantHitData result;
	FindMultiplePeaks peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;


    return true;
}

bool searchForSinglePeakDynamicPedestal()
{
	adcWaveform adcData;
	unsigned int adc[8] = {124, 92, 139, 123, 107, 93, 79, 75};
	adcData = adc;

	resultantHitData result;
	FindMultiplePeaks peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;

    return true;
}

bool searchForDoublePeak()
{
	adcWaveform adcData;
	unsigned int adc[8] = {67, 64, 130, 189, 159, 126, 129, 599};
	adcData = adc;

	resultantHitData result;
	FindMultiplePeaks peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;

    return true;
}

bool searchForDoublePeakConstantPedestal()
{
	adcWaveform adcData;
	unsigned int adc[8] = {80, 72, 213, 388, 350, 536, 530, 384};
	adcData = adc;

	resultantHitData result;
	FindMultiplePeaks peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;

    return true;
}

bool searchForDoublePeakDynamicPedestal()
{
	adcWaveform adcData;
	unsigned int adc[8] = {131, 99, 102, 234, 266, 210, 241, 365};
	adcData = adc;

	resultantHitData result;
	FindMultiplePeaks peak(initParams);
	peak.process(adcData, result);

    std::cout << "result[0]._peakTime" << result[0]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[0]._peakHeight << std::endl;
	std::cout << "result[0]._peakTime" << result[1]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[1]._peakHeight << std::endl;
 	std::cout << "result[0]._peakTime" << result[2]._peakTime << std::endl;
    std::cout << "result[0]._peakHeight" << result[2]._peakHeight << std::endl;

    return true;
}
**/

