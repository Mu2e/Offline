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
		TString prDataFileName = "protonData" + shapingTime + "ns_" + numberOfSample + "DoublePeak6.root";	
		getElectronData(oldtree,shapingTime, elDataFileName);
		getProtonData(oldtree,shapingTime,prDataFileName);
		TString newElFileName = "electronFit" + shapingTime + "ns_" + numberOfSample + "DoublePeak6.root";
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
	

	std::vector<unsigned int> *adc = new std::vector<unsigned int>;
	float mcenergy;
	float mctrigenergy;

	Double_t qMeasurementTimes[numberOfSamples];
	Double_t qAdc[numberOfSamples];
	TGraphErrors graph;
	int qFuncNum;
	float qMcenergy;
	float qMctrigenergy;
	double func7param[5];
	double func8param[6];


	treeData->SetBranchAddress("adc",&adc);
	treeData->SetBranchAddress("mcenergy",&mcenergy);
	treeData->SetBranchAddress("mctrigenergy",&mctrigenergy);
	TFile *newfile = new TFile(newFileName,"RECREATE");

	TTree * convolvedFitTree = new TTree("convolvedFitTree","convolvedFitTree");

	convolvedFitTree->Branch("qAdc",qAdc,"qAdc[" + numSample + "]/D");
	convolvedFitTree->Branch("qMeasurementTimes",qMeasurementTimes,"qMeasurementTimes[" + numSample + "]/D");
	convolvedFitTree->Branch("graph",&graph);
	convolvedFitTree->Branch("func7param",func7param,"func7param[5]/D");
	convolvedFitTree->Branch("func8param",func8param,"func8param[6]/D");
	convolvedFitTree->Branch("qMcenergy",&qMcenergy);
	convolvedFitTree->Branch("qFuncNum",&qFuncNum);
	convolvedFitTree->Branch("qMctrigenergy",&qMctrigenergy);

	for (int trial = 0; trial < 10000; ++trial) //treeData->GetEntries(); ++trial)
	{
		if (trial % 400 == 0)
		{std::cout << trial << endl;}
		treeData->GetEntry(trial);

		// convert adc data from int to double
		for (int k = 0; k < numberOfSamples; ++k){qAdc[k] = (Double_t) (*adc)[k] - 64.0; qMeasurementTimes[k] = (Double_t) 20.0 * k;}

		TF1 *func1 = new TF1("fittingFunction7", fittingFunction7, 0.0, numberOfSamples*20.0, 5);
		const Double_t param1_1 = TMath::Max(TMath::MaxElement(numberOfSamples,qAdc) /0.015, 1000.0);
		const Double_t param1_2 = TMath::Max(qAdc[0],0.0);
		func1->SetParameters(35.0, param1_1, param1_2, 10.0,shapingTime);
		func1->SetParLimits(2,0.0,1000.0);
		func1->SetParLimits(0,0.0,100.0);
		func1->SetParLimits(1,1000.0, 1.0e9);
		func1->SetParLimits(3,0.0,30.0);
		func1->FixParameter(4,shapingTime);

		TF1 *func2 = new TF1("fittingFunction8", fittingFunction8, 0, numberOfSamples*20.0, 6);

		const Double_t param2_1 = TMath::Max(TMath::MaxElement(numberOfSamples,qAdc)/0.015, 1000.0);
		const Double_t param2_2 = qAdc[0];
		const Double_t param2_4 = TMath::Max(TMath::MaxElement(numberOfSamples,qAdc)/0.015, 1000.0);

		func2->SetParameters(35.0, param2_1, param2_2, 50.0 ,param2_4, shapingTime);
		func2->SetParLimits(0,20.0,100.0);
		func2->SetParLimits(1,1000.0,1.0e9);
		func2->SetParLimits(3,25.0,80.0);
		func2->SetParLimits(4,1000.0,1.0e9);
		func2->FixParameter(5,shapingTime); 

		qMcenergy = mcenergy;
		qMctrigenergy = mctrigenergy;


		Double_t errorY[numberOfSamples]; 
		Double_t errorX[numberOfSamples];
		for (int i = 0; i < numberOfSamples; ++i){ errorY[i] = 3.0; errorX[i] = 0.0;}

		graph = TGraphErrors(numberOfSamples,qMeasurementTimes,qAdc,errorX,errorY);

		graph.Fit(func2,"QN");
	    graph.Fit(func1,"QN"); 

		if ( (func2->GetChisquare()*3.0) < func1->GetChisquare())
			{
				qFuncNum = 8;}
		else {
			qFuncNum = 7;
		}

		for (int i = 0; i < 5; ++i){
			func7param[i] = func1->GetParameter(i);
		}
		for (int i = 0; i < 6; ++i){
			func8param[i] = func2->GetParameter(i);
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


