{
	TFile f("electronFit25ns_8DoublePeak6.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");

	float mctrigenergy;
	double qAdc[8];
	double func7param[5];

	treeData->SetBranchAddress("func7param",func7param);
	treeData->SetBranchAddress("qAdc",qAdc);
	treeData->SetBranchAddress("qMctrigenergy",&mctrigenergy);

	Double_t maxElem[treeData->GetEntries()];
	Double_t scaling[treeData->GetEntries()];

	for (int i = 0; i < treeData->GetEntries(); ++i)
	{
		treeData->GetEntry(i);
		maxElem[i] = TMath::MaxElement(8,qAdc);
		scaling[i] = func7param[1];
	}

	TGraph *gr = new TGraph(treeData->GetEntries(),scaling,maxElem);
}