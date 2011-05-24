//
// header for Analyzer.C

// $Id: Analyzer.h,v 1.9 2011/05/24 17:19:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:19:03 $
//
// Original author KLG
//

#include <TROOT.h>
#include <TString.h>
#include <TH1F.h>
#include <TH1.h>
#include <TNtuple.h>
#include <TCanvas.h>

#include <vector>

// avoid "complicated" types in the interface
// one may want to revisit types relying on boost libraries

/* #include "MCDataProducts/inc/SimParticleCollection.hh" */
/* #include "MCDataProducts/inc/GenParticleCollection.hh" */
/* #include "MCDataProducts/inc/StepPointMCCollection.hh" */
/* #include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh" */

/* #include "art/Persistency/Common/Wrapper.h" */
/* #include "art/Framework/Core/Event.h" */

class Analyzer {

 public:

  typedef mu2e::SimParticleCollection::key_type key_type;

//   cint may not be able to handle more complicated defaults in the Dict ...
//   we can't use  TString file="data_03.root",  so we use char const *

  Analyzer (char const * file,
            ULong64_t maxevent = 1000000,
            ULong64_t maxFullPrint = 2,
            Double_t minEnergy = 0.001,
            char const * cformat = "png"
            );

  //  ~Analyzer();

  void begin();
  void analyze();
  void plot();
  void write();

  void printOutCanvases();

  TCanvas* prepareNextCanvas( Int_t nx = 1, Int_t ny = 1,
                              Int_t const logx = 0, Int_t const logy = 0,
                              Int_t const gridx = 1, Int_t const gridy = 1);
  void Analyzer::plotHist(TH1* hist, char const * opt="");
  void Analyzer::plotNHist(std::vector<TH1*>, char const * opt="");
  void Analyzer::plotNT(const char* nts, char const * cut="", char const * opt="");

/* void Analyzer::doLTracker(art::EventAuxiliary*                              EventAuxiliaryWrppd, */
/*                           art::Wrapper<mu2e::StepPointMCCollection>*        StepPointMCWrppd, */
/*                           art::Wrapper<mu2e::GenParticleCollection>*        GenParticleWrppd, */
/*                           art::Wrapper<mu2e::SimParticleCollection>*        SimParticleWrppd, */
/*                           art::Wrapper<mu2e::PhysicalVolumeInfoCollection>* PhysicalVolumeInfoWrppd); */

 private:

  static Long_t const canvasOriginX = 10;
  static Long_t const canvasOriginY = 10;

  static Long_t const canvasWX = 600;
  static Long_t const canvasWY = 600;

  static Long_t const canvasSpace  =  35;
  static Long_t const canvasShiftX =  25;
  static Long_t const canvasShiftY =  25;

  TString _file;
  ULong64_t _mevent;

  // Module label of the g4 module that made the hits.
  //  TString _g4ModuleLabel;

  // Cut on the minimum energy.
  Double_t _minimumEnergy;

  // Limit on number of events for which there will be full printout.
  ULong64_t _maxFullPrint;

  // Number of events analyzed.
  ULong64_t _nAnalyzed;

  // Name for output canvases/files
  TString _outputFileNamePrefix;

  // file format of root canvases
  TString _canvasPrintFormat;

 public:

  // Pointers to histograms to be filled.
  TH1F* _hRadius;
  TH1F* _hEnergyDep;
  TH1F* _hTime;
  TH1F* _hMultiplicity;
  TH1F* _hDriftDist;
  TH1F* _hxHit;
  TH1F* _hyHit;
  TH1F* _hzHit;
  TH1F* _hHitNeighbours;
  TH1F* _hCheckPointRadius;
  TH1F* _hMomentumG4;
  TH1F* _hStepLength;

  TNtuple* _ntup;

  std::vector<TCanvas*> * _canvases;

};
#ifdef __CINT__
#pragma link C++ class Analyzer+;
#endif
