// use 

namespace {
  const char* fn_trkpatrec = "/mu2e/data/users/murat/datasets/tmva_training/e11s5731.tmva_training_trkpatrec.root";
  const char* fn_calpatrec = "/mu2e/data/users/murat/datasets/tmva_training/e11s5731.tmva_training_calpatrec.root";
};


//-----------------------------------------------------------------------------
// Algorithm    : "trkpatrec" or "calpatrec"
// TrainingMode : "chi2d" or "logfcons"
// BkgWeight    : 0,1,2,3,4  (+100 if use  Z)
//-----------------------------------------------------------------------------
int train_tmva(const char* Algorithm = "calpatrec", const char* TrainingMode = "chi2d", int BkgWeight) {

  if ( ! gInterpreter->IsLoaded("TrainTrkQualMurat.C")) {
    gInterpreter->LoadMacro(Form("%s/../source/root/tmva/test/TMVAGui.C",gSystem->Getenv("ROOTSYS")));
    gInterpreter->LoadMacro("CalPatRec/test/TrainTrkQualMurat.C");
  }

  TFile* f;

  TString tmvaName, alg;

  alg = Algorithm;
  
  if      (alg == "trkpatrec") f = TFile::Open(fn_trkpatrec);
  else if (alg == "calpatrec") f = TFile::Open(fn_calpatrec);

  tmvaName  = alg+"_";
  tmvaName += TrainingMode;
  tmvaName += "_";
  tmvaName += "tmva";

  TTree* t = (TTree*) f->Get("tmva_training_tree");

  TrainTrkQualMurat(t,BkgWeight,TrainingMode,tmvaName.Data());
}
