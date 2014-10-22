{
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine(".L dataAnalysis.C+");
	int entry = 0;
	TFile f("electronFit25ns_8DoublePeak6.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");
	TGraphErrors *gr;
	TF1 *func;
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

	std::vector<Float_t> tPeak;
	std::vector<Float_t> adcPeak;
	findPeaks(gr,tPeak,adcPeak,2.0);


}