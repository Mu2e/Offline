//Exampel macro to read scoring plane data in art file
//
//See ReadConfig function to change mesh name and art file name
//

#include <TH2F.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TKey.h>
#include <TROOT.h>
#include <Math/Vector3D.h>






// Internal variables
int nBinsX(0);
int nBinsY(0);
int nBinsZ(0);
ROOT::Math::XYZVectorF halfSize(0,0,0);
ROOT::Math::XYZVectorF center(0,0,0);

namespace mu2e{

  struct ScorerConfigSummary{
    std::string name_;
    unsigned nbinsX_;
    unsigned nbinsY_;
    unsigned nbinsZ_;
    ROOT::Math::XYZVectorF halfSize_;
    ROOT::Math::XYZVectorF center_;
  };

  struct ScorerSummary{
    Int_t ix_;
    Int_t iy_;
    Int_t iz_;
    Int_t entries_;
    Float_t total_;
    Float_t totalSqr_;
  };

  struct GenEventCount{
    unsigned long count_;
  };
}


// ---------------------------------------------------------------------
void readMeshConfig(std::string filename, std::string meshName){

   TTreeReader fReader;

   TTreeReaderValue<mu2e::GenEventCount> eventCounter        = {fReader, "mu2e::GenEventCount_genCounter__POTnocuts.obj"};
   TTreeReaderArray<mu2e::ScorerConfigSummary> configSummary = {fReader, "mu2e::ScorerConfigSummarys_g4run__POTnocuts.obj"};

   std::string name = "mu2e::ScorerSummarys_g4run_"+meshName+"TrackCounter_POTnocuts.obj";
   TTreeReaderArray<mu2e::ScorerSummary> scoresCounter  = {fReader, name.c_str()};
   name = "mu2e::ScorerSummarys_g4run_"+meshName+"DoseDeposit_POTnocuts.obj";
   TTreeReaderArray<mu2e::ScorerSummary> scoresDoseDep  = {fReader, name.c_str()};
   name = "mu2e::ScorerSummarys_g4run_"+meshName+"CellFlux_POTnocuts.obj";
   TTreeReaderArray<mu2e::ScorerSummary> scoresCellFlux = {fReader, name.c_str()};
   name = "mu2e::ScorerSummarys_g4run_"+meshName+"DoseEffective_POTnocuts.obj";
   TTreeReaderArray<mu2e::ScorerSummary> scoresDoseEff = {fReader, name.c_str()};
   name = "mu2e::ScorerSummarys_g4run_"+meshName+"DelayedDose_POTnocuts.obj";
   TTreeReaderArray<mu2e::ScorerSummary> scoresDelayed = {fReader, name.c_str()};

   //Need the first file to read the mesh configurations
   TFile f0(filename.c_str());
   TTree* tree0 = (TTree*) f0.Get("SubRuns");
   fReader.SetTree(tree0);

   fReader.SetEntry(0);
   for (const auto& config : configSummary) {
     if (config.name_.compare(meshName) != 0) continue;
     nBinsX   = config.nbinsX_;
     nBinsY   = config.nbinsY_;
     nBinsZ   = config.nbinsZ_;
     halfSize = config.halfSize_/1000;
     center   = config.center_/1000;
   }
   cout<<"Nbins xyz = "<<nBinsX<<","<<nBinsY<<","<<nBinsZ<<endl;
   cout<<"HalfSize (mm) =("<<halfSize.x()*1000<<","<<halfSize.y()*1000<<","<<halfSize.z()*1000<<")"<<endl;
   cout<<"Center (mm) =("<<center.x()*1000<<","<<center.y()*1000<<","<<center.z()*1000<<")"<<endl;

   f0.Close();
}


// ---------------------------------------------------------------------
void FillHisto(const mu2e::ScorerSummary& score, TH2F& histo){

   float x = (2*halfSize.x())/float(nBinsX)*(score.ix_ +0.5)-halfSize.x() + center.x();
   float y = (2*halfSize.y())/float(nBinsY)*(score.iy_ +0.5)-halfSize.y() + center.y();
   float z = (2*halfSize.z())/float(nBinsZ)*(score.iz_ +0.5)-halfSize.z() + center.z();
   float val = score.total_;

   histo.Fill(z,x,val);
}












void ReadScoring(){

   const std::string meshName("WorldMesh");
   const std::string filename("test.art");


   //Read mesh configuration
   readMeshConfig(filename, meshName);
   TH2F histCounterPOT ("histCounterPOT", "Track Counter (N_{trk}/POT);z (m);x (m)",       nBinsZ, -halfSize.z(), halfSize.z(), nBinsX, -halfSize.x(), halfSize.x());

   //Declare reader
   TTreeReader fReader;
   TTreeReaderValue<mu2e::GenEventCount> eventCounter        = {fReader, "mu2e::GenEventCount_genCounter__POTnocuts.obj"};
   TTreeReaderArray<mu2e::ScorerConfigSummary> configSummary = {fReader, "mu2e::ScorerConfigSummarys_g4run__POTnocuts.obj"};

   std::string name = "mu2e::ScorerSummarys_g4run_"+meshName+"TrackCounter_POTnocuts.obj";
   TTreeReaderArray<mu2e::ScorerSummary> scoresCounter  = {fReader, name.c_str()};
   //add other scorers here (or automate this)


   //Loop over entries and process scoring data
   TFile f1(filename.c_str());
   TTree* tree1 = (TTree*) f1.Get("SubRuns");
   if (tree1 == nullptr) return;
   fReader.SetTree(tree1);

   for (int ievent=0;ievent<tree1->GetEntriesFast();++ievent){
      fReader.SetEntry(ievent);
      for (auto& score : scoresCounter)  FillHisto(score, histCounterPOT);
   }


   //Do what you wish with the histogram
}
