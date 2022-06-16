//
// An EDAnalyzer module that reads the Trigger Info
//
// Original author G. Pezzullo
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "canvas/Utilities/InputTag.h"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// Conditions
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"

// MC dataproducts
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

// ROOT
#include "TH1F.h"
#include "TH2F.h"

#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"

#include <cmath>
// #include <iostream>
#include <string>
// #include <map>
#include <vector>

namespace mu2e {

class EvalWeightedEventCounts : public art::EDAnalyzer {

public:
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  struct Config {
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Debug Level"), 0};
    fhicl::Atom<art::InputTag> evtWeightTag{
        Name("protonBunchIntensity"), Comment("Proton Bunch Intenity "), "protonBunchIntensity"};
    fhicl::Atom<float> nProcess{Name("nEventsProcessed"), Comment("Events Processed"), 1.};
    fhicl::Atom<std::string> tagName{Name("tagName"), Comment("Tag Name"), "tagName"};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  enum {
    kNTrigInfo = 40,
    kNTrackTrig = 20,
    kNTrackTrigVar = 30,
    kNHelixTrig = 40,
    kNHelixTrigVar = 130,
    kNCaloCalib = 5,
    kNCaloCalibVar = 30,
    kNCaloOnly = 5,
    kNCaloOnlyVar = 30,
    kNOcc = 100,
    kNOccVar = 100
  };

  struct occupancyHist_ {
    TH1F* _hOccInfo[kNOcc][kNOccVar];
    TH2F* _h2DOccInfo[kNOcc][kNOccVar];

    occupancyHist_() {
      for (int i = 0; i < kNOcc; ++i) {
        for (int j = 0; j < kNOccVar; ++j) {
          _hOccInfo[i][j] = NULL;
          _h2DOccInfo[i][j] = NULL;
        }
      }
    }
  };

  explicit EvalWeightedEventCounts(const Parameters& config);
  virtual ~EvalWeightedEventCounts() {}

  virtual void beginJob();
  virtual void endJob();

  // This is called for each event.
  virtual void analyze(const art::Event& e);
  virtual void beginRun(const art::Run& run);

  void bookHistograms();
  void bookOccupancyInfoHist(art::ServiceHandle<art::TFileService>& Tfs, occupancyHist_& Hist);

private:
  int _diagLevel;
  art::InputTag _evtWeightTag;

  float _nProcess;
  std::string _tagName;
  double _nPOT;

  occupancyHist_ _occupancyHist;

  const mu2e::Tracker* _tracker;

  // the following pointer is needed to navigate the MC truth info of the strawHits
  const art::Event* _event;

  float _minPOT, _maxPOT;
  double _wgSum[3] = {0};
};

EvalWeightedEventCounts::EvalWeightedEventCounts(const Parameters& config) :
    art::EDAnalyzer(config), _diagLevel(config().diagLevel()), _evtWeightTag(config().evtWeightTag()),
    _nProcess(config().nProcess()), _tagName(config().tagName()), _minPOT(1e6), _maxPOT(4e8) {}

void EvalWeightedEventCounts::bookHistograms() {
  art::ServiceHandle<art::TFileService> tfs;

  bookOccupancyInfoHist(tfs, _occupancyHist);
}

void EvalWeightedEventCounts::bookOccupancyInfoHist(art::ServiceHandle<art::TFileService>& Tfs,
                                                    occupancyHist_& Hist) {

  int index_last = 0;
  art::TFileDirectory occInfoDir = Tfs->mkdir("occInfoGeneral");
  Hist._hOccInfo[index_last][0] = occInfoDir.make<TH1F>(
      Form("hInstLum_%i", index_last), "distrbution of instantaneous lum; p/#mu-bunch", 1000,
      _minPOT, _maxPOT);
  Hist._hOccInfo[index_last][1] = occInfoDir.make<TH1F>(
      Form("hInstLum2B_%i", index_last),
      "distrbution of instantaneous lum re-weighted in 2 batch-mode; p/#mu-bunch", 1000, _minPOT,
      _maxPOT);
}

//--------------------------------------------------------------------------------//
void EvalWeightedEventCounts::beginJob() { bookHistograms(); }

//--------------------------------------------------------------------------------//
void EvalWeightedEventCounts::endJob() {

  //    evalTriggerRate();
  printf("[%s::enJob] totWg    = %10.3e\n", _tagName.c_str(), _wgSum[0]);
  printf("[%s::enJob] totEvtB1 = %10.3f\n", _tagName.c_str(), _wgSum[1]);
  printf("[%s::enJob] totEvtB2 = %10.3f\n", _tagName.c_str(), _wgSum[2]);
}

//--------------------------------------------------------------------------------

//================================================================
void EvalWeightedEventCounts::beginRun(const art::Run& run) {
  // get bfield
  // GeomHandle<BFieldManager> bfmgr;
  // GeomHandle<DetectorSystem> det;
  // CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
  // _bz0 = bfmgr->getBField(vpoint_mu2e).z();

  // mu2e::GeomHandle<mu2e::Tracker> th;
  // _tracker  = th.get();
}

//--------------------------------------------------------------------------------
void EvalWeightedEventCounts::analyze(const art::Event& event) {

  // get the number of POT
  double lumi = -1.;
  art::Handle<ProtonBunchIntensity> evtWeightH;
  event.getByLabel(_evtWeightTag, evtWeightH);
  if (evtWeightH.isValid()) {
    lumi = (double)evtWeightH->intensity();
  }

  // 1 batch mode
  _wgSum[0] += lumi;

  const static double mean_b1 = 1.6e7;
  const static double mean_b2 = 3.9e7;
  const static double sigma = 0.7147;
  const static double mub1 = log(mean_b1) - 0.5 * sigma * sigma;
  const static double mub2 = log(mean_b2) - 0.5 * sigma * sigma;
  const static double xlognorm_norm_b1 =
      6.293492e-8; // evaluated in ROOT for cut_off = 1.2e8,
                   // above mub(1/2) and sigma for x*log-normal(x)
  const static double xlognorm_norm_b2 = 2.887949e-8; // evaluated in ROOT
  const static double cut_off_norm_b1 =
      ROOT::Math::lognormal_cdf(1.2e8, mub1, sigma); // Due to max cutoff in generation
  const static double cut_off_norm_b2 =
      ROOT::Math::lognormal_cdf(1.2e8, mub2, sigma); // Due to max cutoff in generation
  // if(mode == 1) {
  double p1 = lumi * xlognorm_norm_b1 * cut_off_norm_b1;
  _wgSum[1] += p1;
  //   return p1;
  // }
  _occupancyHist._hOccInfo[0][0]->Fill(lumi, p1);

  double p2 = lumi * xlognorm_norm_b2 * cut_off_norm_b2;
  double oneBatchWg = ROOT::Math::lognormal_pdf(lumi, mub1, sigma) / cut_off_norm_b1;
  double twoBatchWg = ROOT::Math::lognormal_pdf(lumi, mub2, sigma) / cut_off_norm_b2;

  p2 *= twoBatchWg / oneBatchWg; // sample generated using one batch mode

  _occupancyHist._hOccInfo[0][1]->Fill(lumi, p2);
  _wgSum[2] += p2;
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EvalWeightedEventCounts);
