#include <sys/types.h>
#include <sys/stat.h>


void TMVACaloNN(const char* files="bkgfiles.txt", std::string basedir="caloNN")
{

    TChain* mytree = new TChain("CaloNNDiag/Calo");

    std::ifstream fs(files,std::ifstream::in);
    char file[100];
    fs.getline(file,100);
    while(fs.good()){
      std::cout << "Adding file " << file << " to training sample" << std::endl;
      mytree->Add(file);
      fs.getline(file,100);
    }

    // Selection cuts
    TCut sigCut("cluConv==1"); 
    TCut bkgCut("cluConv==0"); 

    // Variables to be used for the classifier
    std::vector<std::string> varnames, vardescrip ;
    varnames.push_back("cluEnergy");
    varnames.push_back("cluCogR");
    varnames.push_back("cluNcrys");
    varnames.push_back("cluE1");
    varnames.push_back("cluE2");
    varnames.push_back("cluE9");
    varnames.push_back("cluE25");
    varnames.push_back("cluEout");
  //  varnames.push_back("cluEin");
    vardescrip.push_back("Cluster energy");
    vardescrip.push_back("Cluster radius"); 
    vardescrip.push_back("Number of crystal");
    vardescrip.push_back("Energy Central crystal");
    vardescrip.push_back("Energy E1+E2");
    vardescrip.push_back("Energy E9");
    vardescrip.push_back("Energy E25");
    vardescrip.push_back("Energy out");
 //   vardescrip.push_back("Energy in");
    
    if (varnames.size() != vardescrip.size()) {cout<<"Variable name and description length mismatch, exit"; return;} 


    int dstat = mkdir(basedir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (!dstat) std::cout << "Directory " << basedir << " created " << std::endl;
    std::string rootfile = basedir + std::string("/TMVAoutput.root");
    TFile* outputFile = TFile::Open( rootfile.c_str(), "RECREATE" );


    TMVA::DataLoader *dataloader = new TMVA::DataLoader(basedir);
    for(unsigned ivar=0;ivar<varnames.size();++ivar){
      std::cout << "Adding variable " << varnames[ivar]<< " description: " << vardescrip[ivar] << std::endl;
      dataloader->AddVariable(varnames[ivar],vardescrip[ivar],'F');
    }
    dataloader->AddSignalTree    ( mytree, 1.0);
    dataloader->AddBackgroundTree( mytree, 1.0);
    dataloader->PrepareTrainingAndTestTree(sigCut, bkgCut,"SplitMode=Random:!V:SplitSeed=5309" ); 


    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPRelu", "H:!V:NeuronType=ReLU:VarTransform=N:HiddenLayers=N" );
    //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();    
 
    outputFile->Close();
    delete factory;
}

