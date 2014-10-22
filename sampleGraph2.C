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
	TFile f("electronFit25ns_8DoublePeak6Uniform.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");
	TGraphErrors *gr = new TGraphErrors();
	float mcenergy;
	float mctrigenergy;
	int option;
	TF1 *func = new TF1();
	treeData->SetBranchAddress("graph",&gr);
	treeData->SetBranchAddress("fittingFunction",&func);
	treeData->SetBranchAddress("qMcenergy",&mcenergy);
	treeData->SetBranchAddress("qMctrigenergy",&mctrigenergy);
	treeData->SetBranchAddress("option",&option);
	treeData->GetEntry(entry);

	TFile g("electronFit25ns_8DoublePeak6Gaussian.root");
	TTree *treeDataGauss = (TTree*) gDirectory->Get("convolvedFitTree");


	TCanvas* c1 = new TCanvas("c1", "ElectronFits", 400, 300);
	c1->Divide(2,1);

	const double numberOfSamples = 8.0;
	TF1 *fittingFunction = new TF1();
	TF1 *fittingFunctionGauss = new TF1();

	TF1 *funcGauss = new TF1();
	TGraphErrors *grGauss = new TGraphErrors();
	treeDataGauss->SetBranchAddress("fittingFunction",&funcGauss);
	treeDataGauss->SetBranchAddress("graph",&grGauss);
	treeDataGauss->GetEntry(entry);

	if (option == 1)
	{
		fittingFunction = new TF1("fittingFunction", dynamicPedestal, 0.0, numberOfSamples*20.0, 2);
		fittingFunctionGauss = new TF1("fittingFunctionGauss", dynamicPedestal, 0.0, numberOfSamples*20.0, 2);
	}
	if (option == 2)
	{
		fittingFunction = new TF1("fittingFunction",fittingFunction7Uniform,0.0,numberOfSamples*20.0,5);
		fittingFunctionGauss = new TF1("fittingFunctionGauss",fittingFunction7,0.0,numberOfSamples*20.0,5);
	}
	if (option == 3)
	{
		fittingFunction = new TF1("fittingFunction", fittingFunction8,0.0,numberOfSamples*20.0, 6);
		fittingFunctionGauss = new TF1("fittingFunctionGauss", fittingFunction8,0.0,numberOfSamples*20.0, 6);
	}
	if (option == 4)
	{
		fittingFunction = new TF1("fittingFunction",fittingFunction4Uniform,0.0,numberOfSamples*20.0,6);
		fittingFunctionGauss = new TF1("fittingFunctionGauss",fittingFunction4,0.0,numberOfSamples*20.0,6);
	}
	if (option == 5)
	{
		fittingFunction = new TF1("fittingFunction",fittingFunction10,0.0,numberOfSamples*20.0,6);
		fittingFunctionGauss = new TF1("fittingFunctionGauss",fittingFunction10,0.0,numberOfSamples*20.0,6);
	}

	cout << "option : " << option << endl;
	cout << "par[0]" << func->GetParameters()[0] << endl;
	cout << "par[1]" << func->GetParameters()[1] << endl;
	cout << "par[2]" << func->GetParameters()[2] << endl;
	cout << "par[3]" << func->GetParameters()[3] << endl;
	cout << "par[4]" << func->GetParameters()[4] << endl;
	cout << "par[5]" << func->GetParameters()[5] << endl;
	fittingFunction->SetParameters(func->GetParameters());
	fittingFunctionGauss->SetParameters(funcGauss->GetParameters());

	c1->cd(1);
	//fitFunction7->SetParameters(func7param);
	gr->Draw("A*");
	fittingFunction->Draw("same");
	gr->SetTitle("Uniform");

	c1->cd(2);
	grGauss->Draw("A*");
	fittingFunctionGauss->Draw("same");
	grGauss->SetTitle("Gauss");


	//cout << gr->GetX()[0] << endl;
}