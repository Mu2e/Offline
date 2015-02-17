// This code retrieves the data from the tree in convolvedFit3.root and returns and fits a linear polynomial to the data using specified conditions
{
	TFile f("electronFit25ns_8DoublePeak6Uniform.root");
	TTree *treeData = (TTree*) gDirectory->Get("convolvedFitTree");


Float_t qMctrigenergy;
TF1 *func;
treeData->SetBranchAddress("qMctrigenergy",&qMctrigenergy);
treeData->SetBranchAddress("fittingFunction",&func);
//treeData->SetBranchAddress("qAdc",qAdc);

// Use this line to find out how many elements have mcenergy<0.01
//treeData->Draw("qMcenergy","qMcenergy<0.01"); 27652



Double_t scalingFactors[10000];
Double_t mctrigenergies[10000];

for (int i = 0; i < treeData->GetEntries(); ++i)
{
	treeData->GetEntry(i);

	scalingFactors[i] = func->GetParameter(1);
	mctrigenergies[i] = qMctrigenergy;
}


//Fill Tgraph
TGraph *gr = new TGraph(10000,mctrigenergies,scalingFactors);
TF1 *fitFunc = new TF1("fitFunc","[0]*x",0.0,0.01);

gr->Fit(fitFunc);
gr->Draw("A*");

}



// Results

//Chi2                      =  2.74148e+11
//NDf                       =        27650
//p0                        =      791.974   +/-   28.082      
//p1                        =  6.40668e+06   +/-   11004.9  


