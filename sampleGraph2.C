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
	TFile f("electronFit25ns_8DoublePeak6.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");
	TGraphErrors *gr = new TGraphErrors();
	int qFuncNum;
	float mcenergy;
	float mctrigenergy;
	double func7param[5];
	double func8param[6];
	treeData->SetBranchAddress("graph",&gr);
	treeData->SetBranchAddress("func7param",func7param);
	treeData->SetBranchAddress("func8param",func8param);
	treeData->SetBranchAddress("qFuncNum",&qFuncNum);
	treeData->SetBranchAddress("qMcenergy",&mcenergy);
	treeData->SetBranchAddress("qMctrigenergy",&mctrigenergy);
	treeData->GetEntry(entry);

	TCanvas* c1 = new TCanvas("c1", "ElectronFits", 400, 300);
	c1->Divide(2,1);
	c1->cd(1);

	TF1 *fitFunction7 = new TF1("fittingFunction7", fittingFunction7, 0.0, 8.0 * 20.0, 5);
	fitFunction7->SetParameters(func7param);
	gr->Draw("A*");
	fitFunction7->Draw("same");

	c1->cd(2);
	TF1 *fitFunction8 = new TF1("fittingFunction8", fittingFunction8, 0, 8.0 * 20.0, 6);
	fitFunction8->SetParameters(func8param);
	gr->Draw("A*");
	fitFunction8->Draw("same"); 
}