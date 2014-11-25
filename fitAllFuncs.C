#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <sstream>
#include <vector>
#include "dataAnalysis.C"
#include <iostream>

void getElectronData(TTree *oldtree, TString shapingTime, TString dataFileName);
void getProtonData(TTree *oldtree, TString shapingTime, TString dataFileName);
int fitElectronData(TString dataName, TString newFileName, double shapingTime, const int numberOfSamples);


void fitAllFuncs()
{
	const Double_t shapingTimes[1] = {25.0};
	const int numberOfSamples[1] = {8};

	for (int i = 0; i < 1; ++i)
	{
		for (int j = 0; j < 1; ++j)
		{
		TString shapingTime;
		convert2StringInt(shapingTime, shapingTimes[i]);
		TString numberOfSample;
		convert2StringInt(numberOfSample,numberOfSamples[j]);


		// Get File
		TString fileName = "NoXtalk.root";
		TFile f(fileName);
		gDirectory->cd("makeSD");
		TTree *oldtree = (TTree*) gDirectory->Get("sddiag");
		TString elDataFileName = "electronData" + shapingTime + "ns_" + numberOfSample + ".root";
		TString prDataFileName = "protonData" + shapingTime + "ns_" + numberOfSample + ".root";	
		getElectronData(oldtree,shapingTime, elDataFileName);
		getProtonData(oldtree,shapingTime,prDataFileName);
		TString newElFileName = "electronFit" + shapingTime + "ns_" + numberOfSample + "UniformSecond.root";
		cout << shapingTimes[i] << " : " << numberOfSamples[j] << endl;
		fitElectronData(elDataFileName,newElFileName,shapingTimes[i],numberOfSamples[j]);
		TString newPrFileName = "protonFit" + shapingTime + "ns_" + numberOfSample + "UniformSecond.root";
		fitElectronData(prDataFileName,newPrFileName,shapingTimes[i],numberOfSamples[j]);

	}

	}
}


int fitElectronData(TString dataName, TString newFileName, double shapingTime, const int numberOfSamples)
{
	TFile f(dataName);
	cout << dataName << endl;
	TTree *treeData = (TTree*) gDirectory->Get("sddiag"); 

	TString numSample;
	convert2StringInt(numSample, numberOfSamples);

	int option = 0;
	int suboption = 0;
	TF1 *fittingFunction = new TF1();
	TF1 *fittingFunctionTemp = new TF1(); // Needed to temporarily store function data in high chisquare case

	std::vector<unsigned int> *adc = new std::vector<unsigned int>;
	float mcenergy;
	float mctrigenergy;

	Double_t qMeasurementTimes[numberOfSamples];
	Double_t qAdc[numberOfSamples];
	TGraphErrors *graph = new TGraphErrors();
	float qMcenergy;
	float qMctrigenergy;


	treeData->SetBranchAddress("adc",&adc);
	treeData->SetBranchAddress("mcenergy",&mcenergy);
	treeData->SetBranchAddress("mctrigenergy",&mctrigenergy);
	TFile *newfile = new TFile(newFileName,"RECREATE");

	TTree * convolvedFitTree = new TTree("convolvedFitTree","convolvedFitTree");

	convolvedFitTree->Branch("suboption",&suboption);
	convolvedFitTree->Branch("option",&option);
	convolvedFitTree->Branch("fittingFunction",&fittingFunction);
	convolvedFitTree->Branch("qAdc",qAdc,"qAdc[" + numSample + "]/D");
	convolvedFitTree->Branch("qMeasurementTimes",qMeasurementTimes,"qMeasurementTimes[" + numSample + "]/D");
	convolvedFitTree->Branch("graph",&graph);
	convolvedFitTree->Branch("qMcenergy",&qMcenergy);
	convolvedFitTree->Branch("qMctrigenergy",&qMctrigenergy);

	for (int trial = 0; trial < 1e5; ++trial) //treeData->GetEntries(); ++trial)
	{
		suboption = 0;
		if (trial % 1000 == 0)
		{std::cout << trial << endl;}
		treeData->GetEntry(trial);

		qMcenergy = mcenergy;
		qMctrigenergy = mctrigenergy;

		// convert adc data from int to double
		for (int k = 0; k < numberOfSamples; ++k){
			qAdc[k] = (Double_t) (*adc)[k] - 64.0; 
			qMeasurementTimes[k] = (Double_t) 20.0 * k;}

		Double_t errorY[numberOfSamples]; 
		Double_t errorX[numberOfSamples];
		for (int i = 0; i < numberOfSamples; ++i){ errorY[i] = 3.0; errorX[i] = 0.0;}
		graph = new TGraphErrors(numberOfSamples,qMeasurementTimes,qAdc,errorX,errorY);


		std::vector<Float_t> tPeak;
		std::vector<Float_t> adcPeak;
		findPeaks(graph,tPeak,adcPeak,2.0);
			if (tPeak.size() == 0){
			cout << graph->GetX()[0] << endl;
			cout << "fail : " << trial << endl;}

		fittingFunction = new TF1();
		// If we have a dynamic pedestal
		if (tPeak[0] == 0.0)
		{
			if (tPeak.size() == 1)
			{
				option = 1;

				const double Q = TMath::Max(qAdc[0],0.0);
				Double_t newData[8];

				// Subtract dynamic pedestal from data
				for (int i = 0; i < numberOfSamples; ++i)
				{
					newData[i] = qAdc[i] - Q*exp(-qMeasurementTimes[i]/shapingTime);
				}

				const Double_t newAdcPeak = TMath::MaxElement(numberOfSamples,newData);
				const Double_t newTimePeak = TMath::LocMax(numberOfSamples,newData);

				fittingFunction = new TF1("fittingFunction",fittingFunction7Uniformfixed,0.0,numberOfSamples*20.0,5);
				const double timeShift = newTimePeak - shapingTime - 5.0;
				const double scalingFactor = TMath::Max(newAdcPeak /0.015, 1000.0);
				const double sigma = 10.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,Q,sigma,shapingTime);
				fittingFunction->SetParLimits(2,0.0,1000.0);
				fittingFunction->SetParLimits(0,0.0-shapingTime,140.0-shapingTime); //This is a concession at 0.0
				fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				fittingFunction->SetParLimits(3,0.0,50.0);
				fittingFunction->FixParameter(4,shapingTime);
			}
			if (tPeak.size() == 2)
			{
				option = 2;
				fittingFunction = new TF1("fittingFunction",fittingFunction7Uniformfixed,0.0,numberOfSamples*20.0,5);
				const double timeShift = tPeak[1] - shapingTime - 5.0;
				const double scalingFactor = TMath::Max(adcPeak[1] /0.015, 1000.0);
				const double Q = TMath::Max(qAdc[0],0.0);
				const double sigma = 10.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,Q,sigma,shapingTime);
				fittingFunction->SetParLimits(2,0.0,1000.0);
				fittingFunction->SetParLimits(0,0.0-shapingTime,140.0-shapingTime); //This is a concession at 0.0
				fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				fittingFunction->SetParLimits(3,0.0,50.0);
				fittingFunction->FixParameter(4,shapingTime);
			}
			if (tPeak.size() == 3)
			{
				option = 3;
				fittingFunction = new TF1("fittingFunction", fittingFunction8fixed,0.0,numberOfSamples*20.0, 6);
				const double timeShift0 = tPeak[1] - shapingTime - 5.0;
				const double scalingFactor0 = TMath::Max(adcPeak[1] /0.015, 1000.0);
				const double Q = TMath::Max(qAdc[0],0.0);
				const double timeShift1 = tPeak[2] - tPeak[1];
				const double scalingFactor1 = TMath::Max(adcPeak[2] /0.015, 1000.0);

				fittingFunction->SetParameters(timeShift0, scalingFactor0, Q, timeShift1, scalingFactor1, shapingTime);
				fittingFunction->SetParLimits(0,0.0-shapingTime,80.0-shapingTime);
				fittingFunction->SetParLimits(1,1000.0,1.0e9);
				fittingFunction->SetParLimits(3,25.0,80.0);
				fittingFunction->SetParLimits(4,1000.0,1.0e9);
				fittingFunction->FixParameter(5,shapingTime); 
			}
			graph->Fit(fittingFunction,"QN");
		}
		else
		{

			const double confidence = 3.0;
			// Check if pedestal is significantly above 0
			const bool nonZeroPedestal = (qAdc[0] + qAdc[1])*0.5 > 4.0 / sqrt(2.0) * 3.0;

			if (tPeak.size() == 1)
			{
				option = 4;
				fittingFunction = new TF1("fittingFunction",fittingFunction4Uniformfixed,0.0,numberOfSamples*20.0,6);
				const double timeShift = tPeak[0] - shapingTime - 5.0;
				const double scalingFactor = TMath::Max(adcPeak[0] /0.015, 1000.0);
				const double verticalShift = 0.0;
				const double sigma = 10.0;
				const double truncationLevel = 1024.0 - 64.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime);
				fittingFunction->SetParLimits(0,20.0-shapingTime,140.0-shapingTime); 
				fittingFunction->SetParLimits(1,1000.0, 1.0e9);

				if (nonZeroPedestal)
				{
					fittingFunction->SetParLimits(2,-20.0,1000.0);
				}
				else
				{
					fittingFunction->FixParameter(2,verticalShift);
				}
				fittingFunction->SetParLimits(3,0.0,50.0);
				fittingFunction->FixParameter(4,1024.0-64.0); 
				fittingFunction->FixParameter(5,shapingTime);
				graph->Fit(fittingFunction,"QN");
				//graph->Fit(fittingFunction,"QNM");


				// If Chisquare value of single peak is poor take difference in fit and data 
				// to determine a possible second maximum. Sinc tPeak.size() now equals to
				// the data will now be altered to option 5.
				if (fittingFunction->GetChisquare()>500.0)
				{

					fittingFunctionTemp = fittingFunction;

					//cout << "par0 : " << fittingFunction->GetParameter(0) << endl;
					//cout << "par1 : " << fittingFunction->GetParameter(1) << endl;
					//cout << "par2 : " << fittingFunction->GetParameter(2) << endl;
					//cout << "par3 : " << fittingFunction->GetParameter(3) << endl;
					//cout << "par4 : " << fittingFunction->GetParameter(4) << endl;
					//cout << "par5 : " << fittingFunction->GetParameter(5) << endl;

					double adcResidual[8];
					for (int j = 0; j < numberOfSamples; ++j)
					{
						adcResidual[j] = qAdc[j] - fittingFunction->Eval(qMeasurementTimes[j]);
					}
					const double secondPeak = TMath::MaxElement(numberOfSamples,adcResidual);
					const double locSecondPeak = 20.0*TMath::LocMax(numberOfSamples,adcResidual);

					if (locSecondPeak != 0.0 && abs(locSecondPeak-tPeak[0]) > 20.0)
					{
					suboption = 1;

					if (locSecondPeak > tPeak[0] )
					{
						adcPeak.push_back(secondPeak);
						tPeak.push_back(locSecondPeak);
					}
					// why doesn't locSecondPeak != 0.0 do what it's supposed to do??
					else if (locSecondPeak < tPeak[0] )
					{
						adcPeak.insert(adcPeak.begin(),secondPeak);
						tPeak.insert(tPeak.begin(),locSecondPeak);
					}
				}
					//cout << trial << " : " << tPeak[1] - tPeak[0] << endl;
				}
			}


			if (tPeak.size() == 2)
			{
				option = 5;
				fittingFunction = new TF1("fittingFunction",fittingFunction10fixed,0.0,numberOfSamples*20.0,6);
				const double timeShift0 = tPeak[0] - shapingTime - 5.0;
				const double scalingFactor0 = TMath::Max(adcPeak[0] /0.015, 1000.0);
				double verticalShift;
				const double timeShift1 = tPeak[1] - tPeak[0];
				const double scalingFactor1 = TMath::Max(adcPeak[1] /0.015, 1000.0);

				if (nonZeroPedestal)
				{
					verticalShift = 0.5 * (qAdc[0] + qAdc[1]);
					fittingFunction->SetParLimits(2,-20.0,1000.0);
				}
				else
				{
					verticalShift = 0.0;
					fittingFunction->FixParameter(2,0.0);
				}

				fittingFunction->SetParameters(timeShift0, scalingFactor0, verticalShift, timeShift1, scalingFactor1, shapingTime);	
				fittingFunction->SetParLimits(0,20.0-shapingTime-5.0,140.0-shapingTime);
				fittingFunction->SetParLimits(1,1000.0,1.0e9); 
				fittingFunction->SetParLimits(3,15.0,100.0); 
				fittingFunction->SetParLimits(4,1000.0,1.0e9);
				fittingFunction->FixParameter(5,shapingTime);
				graph->Fit(fittingFunction,"QN");

				// If the 1st fit was "better" keep the 1st fit
				if (fittingFunction->GetChisquare() > fittingFunctionTemp->GetChisquare())
				{
					fittingFunction = fittingFunctionTemp;
					option = 4;
				}

				//graph->Fit(fittingFunction,"QNM");
			}
			//graph->Fit(fittingFunction,"QNM");
		}

		convolvedFitTree->Fill();
	}
		convolvedFitTree->Write();
		//convolvedFitTree->Print(); 
		newfile->Close(); 
		return 0;
}

void getElectronData(TTree *oldtree, TString shapingTime, TString dataFileName)
{
	std::vector<unsigned int> * adc = 0;
	float mcenergy = 0;
	float xTime0;
	Char_t xtalk;
	int mcpdg, mcproc;
    int nend;
    float mcmom;
    float mctrigenergy;

   oldtree->SetBranchAddress("mcmom",&mcmom);
   oldtree->SetBranchAddress("xtalk",&xtalk);
   oldtree->SetBranchAddress("mcpdg",&mcpdg);
   oldtree->SetBranchAddress("mcproc",&mcproc);
   oldtree->SetBranchAddress("nend",&nend);
   oldtree->SetBranchAddress("adc",&adc);
   oldtree->SetBranchAddress("mcenergy",&mcenergy);
   oldtree->SetBranchAddress("mctrigenergy",&mctrigenergy);

    TFile *newfile = new TFile(dataFileName,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    for (Long64_t i=0;i<oldtree->GetEntries(); i++) 
    {
      oldtree->GetEntry(i);
      if (mcpdg==11&&mcproc==56&&!xtalk&&nend==2&&mcenergy!=0.0&&mcmom>50.0) 
      	{newtree->Fill();}
   }
   //newtree->Print();
   newtree->AutoSave();
   delete newfile;
}

void getProtonData(TTree *oldtree, TString shapingTime, TString dataFileName)
{
	std::vector<unsigned int> * adc = 0;
	float mcenergy = 0;
	float xTime0;
	Char_t xtalk;
	int mcpdg, mcproc;
    int nend;
    float mcmom;
    float mctrigenergy;

   oldtree->SetBranchAddress("mcmom",&mcmom);
   oldtree->SetBranchAddress("xtalk",&xtalk);
   oldtree->SetBranchAddress("mcpdg",&mcpdg);
   oldtree->SetBranchAddress("mcproc",&mcproc);
   oldtree->SetBranchAddress("nend",&nend);
   oldtree->SetBranchAddress("adc",&adc);
   oldtree->SetBranchAddress("mcenergy",&mcenergy);
   oldtree->SetBranchAddress("mctrigenergy",&mctrigenergy);

    TFile *newfile = new TFile(dataFileName,"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    for (Long64_t i=0;i<oldtree->GetEntries(); i++) 
    {
      oldtree->GetEntry(i);
      if (mcpdg==2212&&mcproc==56&&!xtalk&&nend==2&&mcenergy!=0.0) 
      	{newtree->Fill();}
   }
   //newtree->Print();
   newtree->AutoSave();
   delete newfile;
 }


void convert2String(TString &string, double doubleNum)
{
	int num = (int) doubleNum;
	std::ostringstream convert;
	convert << num; 
	string = convert.str();
}


