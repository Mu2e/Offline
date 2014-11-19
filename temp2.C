{
	TFile f("electronFit25ns_8DoublePeak6Uniform.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");
	TFile g("electronFit25ns_8DoublePeak6UniformSecond.root");
	TTree *treeDataSecond = (TTree*) gDirectory->Get("convolvedFitTree");

	float qMctrigenergy;
	treeData->SetBranchAddress("qMctrigenergy",&qMctrigenergy);

	int suboption;
	TF1 *funcSecond = new TF1();
	TF1 *func = new TF1();


	treeData->SetBranchAddress("suboption",&suboption);
	treeData->SetBranchAddress("fittingFunction",&func);
	treeDataSecond->SetBranchAddress("fittingFunction",&funcSecond);


	TTree * tempTree = new TTree("tempTree","tempTree");
	float chi1,chi2,energy1,energy2;
	int entry;
	float mctrigenergy;

	tempTree->Branch("chi1",&chi1);
	tempTree->Branch("chi2",&chi2);
	tempTree->Branch("energy1",&energy1);
	tempTree->Branch("energy2",&energy2);
	tempTree->Branch("entry",&entry);
	tempTree->Branch("mctrigenergy",&mctrigenergy);



	for (int i = 0; i < 1e4; ++i)
	{	
		treeDataSecond->GetEntry(i);
		treeData->GetEntry(i);
		if (suboption == 1)
		{
			entry = i;
			chi1 = func->GetChisquare();
			chi2 = funcSecond->GetChisquare();
			energy1 = func->GetParameter(1) / 6.8978069e6;
			energy2 = funcSecond->GetParameter(1) / 6.8978069e6;
			tempTree->Fill();
		}
	}





}

