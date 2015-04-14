#include "TrkChargeReco/inc/ComboPeakFit.hh"
#include <iostream>
#include "TrkChargeReco/inc/FitModelRoot.hh"
//#include "TrkChargeReco/inc/FindPeakComparison.hh"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"

using namespace mu2e::TrkChargeReco;
void testSampleData()
{
 	// Load Sample Data 
	TFile f("TrkChargeReco/test/electronData25ns_8.root");
	TTree *treeData = (TTree*) gDirectory->Get("sddiag");

	const ConfigStruct initParams;	

	adcWaveform adcData;
	Float_t mcenergy;
	std::vector<unsigned int> *adc = new std::vector<unsigned int>;
	treeData->SetBranchAddress("adc",&adc);
	treeData->SetBranchAddress("mcenergy",&mcenergy);
 	
	// Choose processing object
	ComboPeakFit peak(initParams);
	
	Int_t peakNum;
	std::vector<Double_t> *peakHeight = new std::vector<Double_t>;
 	std::vector<Double_t> *peakTime = new std::vector<Double_t>;
	Double_t qMcenergy;
	Double_t computedEnergy;
	
	// Create ttree for fit data
	TFile newfile = TFile("TrkChargeReco/test/FitDataElectron.root","RECREATE");
	TTree fitData = TTree("FitTree","FitTree");	
	fitData.Branch("peakNum",&peakNum);
	fitData.Branch("peakHeight",&peakHeight);
	fitData.Branch("mcenergy",&qMcenergy);
	fitData.Branch("peakTime",&peakTime);
	fitData.Branch("computedEnergy",&computedEnergy);

	for (unsigned int i = 0; i < 1e5; ++i)
	{
		resultantHitData initialGuess;
		unsigned int adcDataTemp[8];
		treeData->GetEntry(i);
		qMcenergy = mcenergy;
		
		// Convert adc data to proper format
		for (int j = 0; j < initParams._numSamplesPerHit; ++j)
		{adcDataTemp[j] = (*adc)[j];}
		adcData = adcDataTemp;
	
		if (i % 1000 == 0)
		std::cout << "i" << i << std::endl;
		peak.initialPeakGuess(adcData, initialGuess);
		
		// For fits use this line
	//	resultantHitData result(initialGuess.size());

		// For sum method use this line
		resultantHitData result(1);

		peak.process(adcData, initialGuess, result);
		for (unsigned int k = 0; k < result.size(); ++k)
		{
			peakHeight->push_back(result[k]._peakHeight);
			peakTime->push_back(result[k]._peakTime);
		}
		if (result[0]._peakTime == 0.0)
			computedEnergy = result[1]._peakHeight;
		else
			computedEnergy = result[0]._peakHeight;

		// For sum ADC
		//computedEnergy = result[0]._peakHeight;

		peakNum = result.size();
		fitData.Fill();
		peakHeight->clear();
		peakTime->clear();
	}

		fitData.Write();
		newfile.Close(); 
}



/**void sampleGraph(const int num)
{
 	// Load Sample Data 
	TFile f("TrkChargeReco/test/protonData25ns_8.root");
	TTree *treeData = (TTree*) gDirectory->Get("sddiag");

	const ConfigStruct initParams;	

	adcWaveform adcData;
	Float_t mcenergy;
	std::vector<unsigned int> *adc = new std::vector<unsigned int>;
	treeData->SetBranchAddress("adc",&adc);
	treeData->SetBranchAddress("mcenergy",&mcenergy);
	FindMultiplePeaks peak(initParams);

	Int_t peakNum;
	std::vector<Double_t> *peakHeight = new std::vector<Double_t>;
 	std::vector<Double_t> *peakTime = new std::vector<Double_t>;
	Double_t qMcenergy;
	Double_t computedEnergy;


	Double_t ADCdata[8];
	
		resultantHitData initialGuess;
		unsigned int adcDataTemp[8];
		treeData->GetEntry(num);
		qMcenergy = mcenergy;
		for (int j = 0; j < initParams._numSamplesPerHit; ++j)
		{adcDataTemp[j] = (*adc)[j];
		ADCdata[j] = (*adc)[j];}
		adcData = adcDataTemp;
		peak.initialPeakGuess(adcData, initialGuess);
		resultantHitData result(initialGuess.size());
		peak.process(adcData, initialGuess, result);
	//	for (unsigned int k = 0; k < result.size(); ++k)
	//	{
	//		peakHeight->push_back(result[k]._peakHeight);
	//		peakTime->push_back(result[k]._peakTime);
	//	}
		if (result[0]._peakTime == 0.0)
			computedEnergy = result[1]._peakHeight;
		else
			computedEnergy = result[0]._peakHeight;

		peakNum = result.size();


		Double_t measurementTimes[8];
		for (int i = 0; i < 8; ++i)
		{
			measurementTimes[i] = 20.0*i;
		}
		TGraph *graph = new TGraph(8,measurementTimes,ADCdata);
		graph->Draw("A*");
		std::cout << "numPeaks" << result.size() << std::endl;
		std::cout << "peakTime" <<  result[0]._peakTime << std::endl;
		std::cout << "peakHeight" << result[0]._peakHeight << std::endl;
		std::cout << "computed Energy : " << computedEnergy << std::endl;
		std::cout << "initialGuess peakheight : " << initialGuess[0]._peakHeight << std::endl;
		std::cout << "initialGuess peaktime : " << initialGuess[0]._peakTime << std::endl;

		
}**/
