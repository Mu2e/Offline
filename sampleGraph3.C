#include <vector>
#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TRoot.h"
#include "dataAnalysis.C"


void sampleGraph(Long64_t entry)
{
	TFile f("electronFit25ns_8UniformFirst.root");

	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");
	TGraphErrors *gr = new TGraphErrors();
	float mcenergy;
	float mctrigenergy;
	int option;
	int suboption;
	TF1 *func = new TF1();
	treeData->SetBranchAddress("graph",&gr);
	treeData->SetBranchAddress("fittingFunction",&func);
	treeData->SetBranchAddress("qMcenergy",&mcenergy);
	treeData->SetBranchAddress("qMctrigenergy",&mctrigenergy);
	treeData->SetBranchAddress("option",&option);
	treeData->SetBranchAddress("suboption",&suboption);
	treeData->GetEntry(entry);

	TFile g("electronFit25ns_8UniformSecond.root");
	TTree *treeDataGauss = (TTree*) gDirectory->Get("convolvedFitTree");


	TCanvas* c1 = new TCanvas("c1", "ElectronFits", 400, 300);
	c1->Divide(2,1);

	const double numberOfSamples = 8.0;
	TF1 *fittingFunction = new TF1();
	TF1 *fittingFunctionGauss = new TF1();
	int option2;
	int suboption2;

	TF1 *funcGauss = new TF1();
	TGraphErrors *grGauss = new TGraphErrors();
	treeDataGauss->SetBranchAddress("fittingFunction",&funcGauss);
	treeDataGauss->SetBranchAddress("graph",&grGauss);
	treeDataGauss->SetBranchAddress("option",&option2);
	treeDataGauss->SetBranchAddress("suboption",&suboption2);
	treeDataGauss->GetEntry(entry);

	if (option == 1)
	{
		fittingFunction = new TF1("fittingFunction", fittingFunction7Uniformfixed, 0.0, numberOfSamples*20.0, 5);
		fittingFunctionGauss = new TF1("fittingFunctionGauss", fittingFunction7Uniformfixed, 0.0, numberOfSamples*20.0, 5);
	}
	if (option2 == 1)
		fittingFunctionGauss = new TF1("fittingFunctionGauss", fittingFunction7Uniformfixed, 0.0, numberOfSamples*20.0, 5);
	if (option == 2)
		fittingFunction = new TF1("fittingFunction",fittingFunction7Uniformfixed,0.0,numberOfSamples*20.0,5);
	if (option2 == 2)
		fittingFunctionGauss = new TF1("fittingFunctionGauss",fittingFunction7Uniformfixed,0.0,numberOfSamples*20.0,5);
	if (option == 3)
		fittingFunction = new TF1("fittingFunction", fittingFunction8fixed,0.0,numberOfSamples*20.0, 6);
	if (option2 == 3)
		fittingFunctionGauss = new TF1("fittingFunctionGauss", fittingFunction8fixed,0.0,numberOfSamples*20.0, 6);
	if (option == 4)
		fittingFunction = new TF1("fittingFunction",fittingFunction4Uniformfixed,0.0,numberOfSamples*20.0,6);
	if (option2 == 4)
		fittingFunctionGauss = new TF1("fittingFunctionGauss",fittingFunction4Uniformfixed,0.0,numberOfSamples*20.0,6);
	if (option == 5)
		fittingFunction = new TF1("fittingFunction",fittingFunction10fixed,0.0,numberOfSamples*20.0,6);
	if (option2 == 5)
		fittingFunctionGauss = new TF1("fittingFunctionGauss",fittingFunction10fixed,0.0,numberOfSamples*20.0,6);

	cout << "option : " << option2 << endl;
	cout << "par[0] : " << funcGauss->GetParameters()[0] << endl;
	cout << "par[1] : " << funcGauss->GetParameters()[1] << endl;
	cout << "par[2] : " << funcGauss->GetParameters()[2] << endl;
	cout << "par[3] : " << funcGauss->GetParameters()[3] << endl;
	cout << "par[4] : " << funcGauss->GetParameters()[4] << endl;
	cout << "par[5] : " << funcGauss->GetParameters()[5] << endl;
	cout << "MC Trigenergy" << mctrigenergy << endl;
	cout << "suboption : " << suboption2 << endl;  
	fittingFunction->SetParameters(func->GetParameters());
	fittingFunctionGauss->SetParameters(funcGauss->GetParameters());

	c1->cd(1);
	//fitFunction7->SetParameters(func7param);
	gr->Draw("A*");
	fittingFunction->Draw("same");
	gr->SetTitle("Dynamic Truncation");

	c1->cd(2);
	grGauss->Draw("A*");
	fittingFunctionGauss->Draw("same");
	grGauss->SetTitle("Fixed Truncation");

	//cout << gr->GetX()[0] << endl;
}