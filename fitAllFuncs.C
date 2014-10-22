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
		TString fileName = "xtalk_" + shapingTime + "ns_" + numberOfSample + ".root";
		TFile f(fileName);
		gDirectory->cd("makeSD");
		TTree *oldtree = (TTree*) gDirectory->Get("sddiag");
		TString elDataFileName = "electronData" + shapingTime + "ns_" + numberOfSample + "DoublePeak6.root";
		//TString prDataFileName = "protonData" + shapingTime + "ns_" + numberOfSample + "DoublePeak6.root";	
		//getElectronData(oldtree,shapingTime, elDataFileName);
		//getProtonData(oldtree,shapingTime,prDataFileName);
		TString newElFileName = "electronFit" + shapingTime + "ns_" + numberOfSample + "DoublePeak6Uniform.root";
		cout << shapingTimes[i] << " : " << numberOfSamples[j] << endl;
		fitElectronData(elDataFileName,newElFileName,shapingTimes[i],numberOfSamples[j]);
		//TString newPrFileName = "protonFit" + shapingTime + "ns_" + numberOfSample + "DoublePeak6.root";
		//fitElectronData(prDataFileName,newPrFileName,shapingTimes[i],numberOfSamples[j]);
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
	TF1 *fittingFunction = new TF1();
	

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


	convolvedFitTree->Branch("option",&option);
	convolvedFitTree->Branch("fittingFunction",&fittingFunction);
	convolvedFitTree->Branch("qAdc",qAdc,"qAdc[" + numSample + "]/D");
	convolvedFitTree->Branch("qMeasurementTimes",qMeasurementTimes,"qMeasurementTimes[" + numSample + "]/D");
	convolvedFitTree->Branch("graph",&graph);
	convolvedFitTree->Branch("qMcenergy",&qMcenergy);
	convolvedFitTree->Branch("qMctrigenergy",&qMctrigenergy);

	for (int trial = 0; trial < 10000; ++trial) //treeData->GetEntries(); ++trial)
	{
		if (trial % 400 == 0)
		{std::cout << trial << endl;}
		treeData->GetEntry(trial);

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
				fittingFunction = new TF1("fittingFunction", dynamicPedestal, 0.0, numberOfSamples*20.0, 2);
				const double Q = qAdc[0];
				fittingFunction->SetParameters(Q,shapingTime);
				fittingFunction->SetParLimits(0,0.0,1000.0);
				fittingFunction->FixParameter(1,shapingTime);
			}
			if (tPeak.size() == 2)
			{
				option = 2;
				fittingFunction = new TF1("fittingFunction",fittingFunction7Uniform,0.0,numberOfSamples*20.0,5);
				const double timeShift = tPeak[1] - shapingTime;
				const double scalingFactor = TMath::Max(adcPeak[1] /0.015, 1000.0);
				const double Q = TMath::Max(qAdc[0],0.0);
				const double sigma = 10.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,Q,sigma,shapingTime);
				fittingFunction->SetParLimits(2,0.0,1000.0);
				fittingFunction->SetParLimits(0,0.0-shapingTime,140.0-shapingTime); //This is a concession at 0.0
				fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				fittingFunction->SetParLimits(3,0.0,30.0);
				fittingFunction->FixParameter(4,shapingTime);
			}
			if (tPeak.size() == 3)
			{
				option = 3;
				fittingFunction = new TF1("fittingFunction", fittingFunction8,0.0,numberOfSamples*20.0, 6);
				const double timeShift0 = tPeak[1] - shapingTime;
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
			if (tPeak.size() == 1)
			{
				option = 4;
				fittingFunction = new TF1("fittingFunction",fittingFunction4Uniform,0.0,numberOfSamples*20.0,6);
				const double timeShift = tPeak[0] - shapingTime;
				const double scalingFactor = TMath::Max(adcPeak[0] /0.015, 1000.0);
				const double verticalShift = 0.5 * (qAdc[0] + qAdc[1]);
				const double sigma = 10.0;
				const double truncationLevel = 1024.0 - 64.0;
				fittingFunction->SetParameters(timeShift,scalingFactor,verticalShift,sigma,truncationLevel,shapingTime);
				fittingFunction->SetParLimits(0,20.0-shapingTime,140.0-shapingTime); 
				fittingFunction->SetParLimits(1,1000.0, 1.0e9);
				fittingFunction->SetParLimits(2,-20.0,1000.0); // This bound currently works but is questionable
				fittingFunction->SetParLimits(3,0.0,30.0);
				fittingFunction->FixParameter(4,1024.0-64.0); 
				fittingFunction->FixParameter(5,shapingTime);
			}
			if (tPeak.size() == 2)
			{
				option = 5;
				fittingFunction = new TF1("fittingFunction",fittingFunction10,0.0,numberOfSamples*20.0,6);
				const double timeShift0 = tPeak[0] - shapingTime;
				const double scalingFactor0 = TMath::Max(adcPeak[0] /0.015, 1000.0);
				const double verticalShift = 0.5 * (qAdc[0] + qAdc[1]);
				const double timeShift1 = tPeak[1] - tPeak[0];
				const double scalingFactor1 = TMath::Max(adcPeak[1] /0.015, 1000.0);
				fittingFunction->SetParameters(timeShift0, scalingFactor0, verticalShift, timeShift1, scalingFactor1, shapingTime);
				fittingFunction->SetParLimits(0,20.0-shapingTime,100.0-shapingTime);
				fittingFunction->SetParLimits(1,1000.0,1.0e9);
				fittingFunction->SetParLimits(2,-20.0,1000.0); // This bound currently works but is questionable
				fittingFunction->SetParLimits(3,25.0,100.0);
				fittingFunction->SetParLimits(4,1000.0,1.0e9);
				fittingFunction->FixParameter(5,shapingTime); 
			}
			graph->Fit(fittingFunction,"QN");
		}

		convolvedFitTree->Fill();
	}
		convolvedFitTree->Write();
		convolvedFitTree->Print(); 
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
   newtree->Print();
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
   newtree->Print();
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


