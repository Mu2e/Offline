//
// An EvalWeightedEventCount module that analyzes Trigger Rates
//
// Original author Tausif Hossain
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

// MC dataproducts
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

// ROOT
#include "TH1F.h"

#include "Math/PdfFuncMathCore.h"
#include "Math/ProbFuncMathCore.h"
#include <cmath>
#include <string>
#include <vector>

namespace mu2e {

class EvalWeightedEventCounts : public art::EDAnalyzer {

public:
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  struct Config {
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagLevel"), 0};
    fhicl::Atom<art::InputTag> evtWeightTag{Name("protonBunchIntensity"),
                                            Comment("protonBunchIntensity")};
    fhicl::Atom<std::string> tagName{Name("tagName"), Comment("tagName")};
  };

  using Parameters = art::EDAnalyzer::Table<Config>;

  enum { kNOcc = 2, kNOccVar = 2 };

  struct occupancyHist_ {
    TH1F* _hOccInfo[kNOcc][kNOccVar];

    occupancyHist_() {
      for (int i = 0; i < kNOcc; ++i) {
        for (int j = 0; j < kNOccVar; ++j) {
          _hOccInfo[i][j] = NULL;
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

  std::string _tagName;
  double _nPOT;

  occupancyHist_ _occupancyHist;

  float _minPOT, _maxPOT;
  std::array<double, 3> _wgSum = {0};
};

EvalWeightedEventCounts::EvalWeightedEventCounts(const Parameters& config) :
    art::EDAnalyzer(config), _diagLevel(config().diagLevel()),
    _evtWeightTag(config().evtWeightTag()), _tagName(config().tagName()), _minPOT(1e6),
    _maxPOT(4e8) {}

void EvalWeightedEventCounts::bookHistograms() {
  art::ServiceHandle<art::TFileService> tfs;

  bookOccupancyInfoHist(tfs, _occupancyHist);
}

void EvalWeightedEventCounts::bookOccupancyInfoHist(art::ServiceHandle<art::TFileService>& Tfs,
                                                    occupancyHist_& Hist) {

  int index_last = 0;
  art::TFileDirectory occInfoDir = Tfs->mkdir("occInfoGeneral");
  Hist._hOccInfo[index_last][0] =
      occInfoDir.make<TH1F>(Form("hInstLum_%i", index_last),
                            "distrbution of instantaneous lum; p/pulse", 1000, _minPOT, _maxPOT);
  Hist._hOccInfo[index_last][1] =
      occInfoDir.make<TH1F>(Form("hInstLum2B_%i", index_last),
                            "distrbution of instantaneous lum re-weighted in 2 batch-mode; p/pulse",
                            1000, _minPOT, _maxPOT);
}

//--------------------------------------------------------------------------------//
void EvalWeightedEventCounts::beginJob() { bookHistograms(); }

//--------------------------------------------------------------------------------//
void EvalWeightedEventCounts::endJob() {

  //    evalTriggerRate();  Tag names must be specified in the fcl
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
  auto const evtWeightH = event.getValidHandle<ProtonBunchIntensity>(_evtWeightTag);
  double lumi = evtWeightH->intensity();

  // 1 batch mode
  _wgSum[0] += lumi;

  const static double mean_b1 = 1.6e7; // Mean Proton Pulse Intensity in 1-Batch Mode
  const static double mean_b2 = 3.9e7; // Mean  Proton Pulse Intensity in 2-Batch Mode
  const static double sigma = 0.7147;  // Standard deviation of the x*log-normal(x) of the Proton
                                       // Pulse Intensity distribution
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

  double p1 = lumi * xlognorm_norm_b1 * cut_off_norm_b1;
  _wgSum[1] += p1;
  _occupancyHist._hOccInfo[0][0]->Fill(lumi, p1);

  double p2 = lumi * xlognorm_norm_b2 * cut_off_norm_b2;
  double oneBatchWg = ROOT::Math::lognormal_pdf(lumi, mub1, sigma) / cut_off_norm_b1;
  double twoBatchWg = ROOT::Math::lognormal_pdf(lumi, mub2, sigma) / cut_off_norm_b2;

  p2 *= twoBatchWg / oneBatchWg; // sample generated using one batch mode

  _occupancyHist._hOccInfo[0][1]->Fill(lumi, p2);
  _wgSum[2] += p2;
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EvalWeightedEventCounts)
