{
	gROOT->ProcessLine("#include <vector>");
	gROOT->ProcessLine(".L dataAnalysis.C+");
	TFile f("electronFit25ns_8DoublePeak6.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");
	TGraphErrors *gr;
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

	for (int i = 0; i < treeData->GetEntries(); ++i)
	{
		std::vector <int> sizes2;
		std::vector<int> sizes3;
		std::vector<Float_t> tPeak3;
		std::vector<Float_t> adcPeak3;
		std::vector<Float_t> tPeak2;
		std::vector<Float_t> adcPeak2;
		treeData->GetEntry(i);
		findPeaks(gr, tPeak3, adcPeak3,3.0);
		findPeaks(gr, tPeak2, adcPeak2,2.0);

		if (tPeak2.size() >= 2 && tPeak3.size() == 1)
		{
			if (tPeak2[0] != 0.0)
			{cout << i << "Peak 2 size : " << tPeak2.size() << endl;}
		}
	}

}