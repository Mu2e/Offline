//Macro to train the NN for TimeAndPhiClusterFinder module. Need to use TPCFNN.C to produce the required data
void TPCFTMVATrain(TString cutTuple="")
{
    TMVA::DataLoader *dataloader = new TMVA::DataLoader("NN");
    TFile* outputFile = TFile::Open("TMVAout.root", "RECREATE" );

    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");

    dataloader->AddVariable( "var2",   "Radius",         "mm",       'F' );
    dataloader->AddVariable( "var3",   "Phi",         "mm",       'F' );
    dataloader->AddVariable( "var4",   "Time",         "mm",     'F' );

    TFile *input = TFile::Open("data.root");
    TTree *sig   = (TTree*)input->Get("Data");
    TTree *bkg   = (TTree*)input->Get("Data");

    dataloader->AddSignalTree(sig);
    dataloader->AddBackgroundTree(bkg);

    dataloader->PrepareTrainingAndTestTree("var1==1", "var1==0","SplitMode=Random:!V:SplitSeed=5309" );
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:VarTransform=Norm:HiddenLayers=N:NCycles=300:NeuronType=ReLU:" );

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();
    delete factory;
}
