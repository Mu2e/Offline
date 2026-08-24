// =============================================================================
// extrapolate the MIP MPV for offline calorimeter calibration with cosmics
// 09-Mar-2026 - S.Salamino & S.Giovannella
//
// Starting code for track selection and path normalization:
// CaloCosmicEnergy module from P.Fedeli & S.Giovannella
//
// =============================================================================

// includes

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ROOT includes
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"

namespace mu2e {

class CaloCosmicEnecalib : public art::EDAnalyzer {

public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic level"), 0};
    fhicl::Atom<bool> useMeV{Name("useMeV"), Comment("Set to true for data in MeV, false for ADC"),
                             false};
    fhicl::Atom<int> CutNCryHit{Name("CutNCryHit"),
                                Comment("Minimum number of crystals in the event"), 3};
    fhicl::Atom<float> CutEnergyDep{Name("CutEnergyDep"), Comment("Minimum energy of the hit"),
                                    250.}; // ADC
    fhicl::Atom<float> LYSOcut{Name("LYSOcut"), Comment("Minimum energy of the hit for LYSO"),
                               500.}; // ADC
    fhicl::Atom<float> CutChi2Norm{Name("CutChi2Norm"),
                                   Comment("Maximum chi2 threshold for linear fits"), 2.5};
    fhicl::Atom<art::InputTag> CaloClusterTag{
        Name("CaloClusterTag"), Comment("Tag for Calorimeter cluster collection"), art::InputTag()};
    fhicl::Atom<art::InputTag> CaloHitTag{
        Name("CaloHitTag"), Comment("Tag for Calorimeter hit collection"), art::InputTag()};
    fhicl::Atom<std::string> OutCalibFile{Name("OutCalibFile"),
                                          Comment("Name for output .dat MIP calibration file"),
                                          "calib_parameters.dat"};
  };

  explicit CaloCosmicEnecalib(const art::EDAnalyzer::Table<Config>& config);
  virtual ~CaloCosmicEnecalib(){};
  void beginJob() override;
  void beginRun(art::Run const& run) override;
  void analyze(art::Event const& event) override;
  void endJob() override;

private:
  art::ServiceHandle<art::TFileService> tfs;

  art::TFileDirectory tfdir = tfs->mkdir("All_tracks");
  art::TFileDirectory tfdirv = tfs->mkdir("Vertical_tracks");
  art::TFileDirectory tfdird = tfs->mkdir("Diagonal_tracks");
  art::TFileDirectory tfdirg = tfs->mkdir("General_tracks");
  art::TFileDirectory tfdir_res = tfs->mkdir("MPV_distributions");

  // calo parameters
  static constexpr int nSiPMs = CaloConst::_nSiPMPerCrystal;
  static constexpr int nCrystals = CaloConst::_nCrystalPerDisk;
  static constexpr int nDisks = CaloConst::_nDisk;
  static constexpr int nROchan = nDisks * nCrystals * nSiPMs;
  const Calorimeter* cal; // For calorimeter geometry
  float cryDim;           // Total crystal dimension
  float MaxDxVertical;    // Max Dx value for vertical tracks

  // fcl parameters
  int _diagLevel;
  bool _useMeV;
  int CutNCryHit;     // >
  float CutEnergyDep; // >
  float LYSOcut;      // >
  float CutChi2Norm;  // <
  art::ProductToken<CaloClusterCollection> _caloClusterToken;
  art::ProductToken<CaloHitCollection> _caloHitToken;
  std::ofstream _outputfile;

  // histogram variables
  TH1F* hSiPM[nROchan] = {nullptr};
  TH1F* hSiPMv[nROchan] = {nullptr};
  TH1F* hSiPMd[nROchan] = {nullptr};
  TH1F* hSiPMg[nROchan] = {nullptr};

  TH1F* MPV = nullptr;
  TH1F* MPVv = nullptr;
  TH1F* MPVd = nullptr;
  TH1F* MPVg = nullptr;
  TH1F* hWidth = nullptr;
  TH1F* hSigma = nullptr;

  int hSiPMbins = 330;        // bins for hSiPM histograms
  int lowStatThreshold = 370; // Number of events threshold to rebin histos with low statistics
  int midStatThreshold = 600; // Number of events threshold to rebin histos with medium statistics
  double hSiPMxmax;
  double hMPVxmax;
  double hWxmax;
  double hSxmax;

  TF1* myPoly = nullptr;
  TF1* mylang = nullptr;
  TF1* mygaus = nullptr;

  // helper variables for Langaus fit
  struct LangausResult { // langaus fit parameters
    float mpv = -99.;
    float mpvErr = -99.;
    float width = -99.;
    float widthErr = -99.;
    float sigma = -99.;
    float sigmaErr = -99.;
    float chi2 = -99.;
    int ndf = 0;
    int nev = 0;
    void reset() { *this = LangausResult(); }
  };

  struct FitFlags { // flags to determine fit quality
    bool badWidthLow = false;
    bool badWidthHigh = false;
    bool badSigmaLow = false;
    bool badSigmaHigh = false;
    bool isBad() const { return badWidthLow || badWidthHigh || badSigmaLow || badSigmaHigh; }
    void fitQuality(float width, float sigma, TH1F* hW, TH1F* hS) {
      if (!hW || !hS || width <= 0 || sigma <= 0)
        return;
      badWidthLow = (width < hW->GetMean() - 3 * hW->GetRMS());
      badWidthHigh = (width > hW->GetMean() + 4 * hW->GetRMS());
      badSigmaLow = (sigma < hS->GetMean() - 3 * hS->GetRMS());
      badSigmaHigh = (sigma > hS->GetMean() + 4 * hS->GetRMS());
    }
  };

  std::vector<LangausResult> sipmParams;
  std::vector<LangausResult> sipmParamsv;
  std::vector<LangausResult> sipmParamsd;
  std::vector<LangausResult> sipmParamsg;

  static double langaus(Double_t* x, Double_t* par);
  float findpath(float m, float q, float x, float y);
  LangausResult fitLangaus(TH1F* histo, int refit = 0, int iter = 0, float oldwidth = 0,
                           float oldsigma = 0);

}; // end class

CaloCosmicEnecalib::CaloCosmicEnecalib(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer(config), _diagLevel(config().diagLevel()), _useMeV(config().useMeV()),
    CutNCryHit(config().CutNCryHit()), CutEnergyDep(config().CutEnergyDep()),
    LYSOcut(config().LYSOcut()), CutChi2Norm(config().CutChi2Norm()),
    _caloClusterToken(consumes<CaloClusterCollection>(config().CaloClusterTag())),
    _caloHitToken(consumes<CaloHitCollection>(config().CaloHitTag())),
    _outputfile(config().OutCalibFile()) // opens file.dat
{
  if (!_outputfile.is_open()) {
    throw cet::exception("CaloCosmicEnecal")
        << "ERROR! Cannot open output file " << config().OutCalibFile() << std::endl;
  }
}

void CaloCosmicEnecalib::beginRun(art::Run const& run) {

  art::ServiceHandle<GeometryService> geom;
  if (!geom->hasElement<Calorimeter>()) {
    throw cet::exception("NO-CALO-GEOM") << "CaloEnergy: Calorimeter geometry not found";
  }

  GeomHandle<Calorimeter> ch;
  cal = ch.get();
  cryDim = cal->G4Info().get<double>("crystalXYLength") +
           2. * cal->G4Info().get<double>("wrapperThickness");
  MaxDxVertical = cryDim * 1.1;
}

void CaloCosmicEnecalib::beginJob() {
  if (_diagLevel > 0)
    std::cout << "CaloCosmicEnecalib: Entering beginJob" << std::endl;

  if (!_useMeV) {
    hSiPMxmax = 2000.;
    hMPVxmax = 1000.;
    hWxmax = 200.;
    hSxmax = 400.;
  } else {
    hSiPMxmax = 100.;
    hMPVxmax = 35.;
    hWxmax = 10.;
    hSxmax = 20.;
  }

  for (int iRou = 0; iRou < nROchan; iRou++) {
    hSiPM[iRou] =
        tfdir.make<TH1F>(Form("Channel_%i", iRou), Form("Energy deposited in channel %i", iRou),
                         hSiPMbins, 0, hSiPMxmax);
    hSiPMv[iRou] = tfdirv.make<TH1F>(Form("Channel_%i_v", iRou),
                                     Form("Energy deposited in channel %i - vertical tracks", iRou),
                                     80, 0, hSiPMxmax);
    hSiPMd[iRou] = tfdird.make<TH1F>(Form("Channel_%i_d", iRou),
                                     Form("Energy deposited in channel %i - diagonal tracks", iRou),
                                     80, 0, hSiPMxmax);
    hSiPMg[iRou] = tfdirg.make<TH1F>(Form("Channel_%i_g", iRou),
                                     Form("Energy deposited in channel %i - general tracks", iRou),
                                     80, 0, hSiPMxmax);
  }

  MPV = tfdir_res.make<TH1F>("MPV_distr", "MPV distribution", 400, 0., hMPVxmax);
  MPVv = tfdir_res.make<TH1F>("MPV_v_distr", "MPV distribution only vertical tracks", 1000, 15.,
                              hMPVxmax);
  MPVd = tfdir_res.make<TH1F>("MPV_d_distr", "MPV distribution only diagonal tracks", 1000, 15.,
                              hMPVxmax);
  MPVg = tfdir_res.make<TH1F>("MPV_g_distr", "MPV distribution only general tracks", 1000, 15.,
                              hMPVxmax);
  hWidth = tfdir_res.make<TH1F>("Width_distr", "Width distribution", 150, 0., hWxmax);
  hSigma = tfdir_res.make<TH1F>("Sigma_distr", "Sigma distribution", 300, 0., hSxmax);

  myPoly = new TF1("myPoly", "pol1", -600., 600.);   // polinomial to fit the track
  mylang = new TF1("mylang", langaus, 0., 2100., 4); // langaus function
  mygaus = new TF1("mygaus", "gaus");

} // end begin job

void CaloCosmicEnecalib::analyze(art::Event const& event) {

  auto const& caloClusters = event.getProduct(_caloClusterToken);
  int nCluster = caloClusters.size();

  // Loop over clusters in the event
  for (int iClu = 0; iClu < nCluster; iClu++) {
    if (_diagLevel > 0)
      std::cout << "iClu: " << iClu << std::endl;

    int nhits = 0;    // number of crystal hits passing the selection
    float Dx = 1400.; // to calculate the spread along x

    // creating arrays to store information for hits passing the selection
    std::array<int, nCrystals> IDs;
    std::array<double, nCrystals> Vmax;
    std::array<float, nCrystals> path;
    std::vector<std::pair<float, float>> PosXY;

    auto const& hitsPtrVector = caloClusters[iClu].caloHitsPtrVector();
    int diskID = caloClusters[iClu].diskID();

    // Loop over crystals in the cluster
    for (auto const& hitsPtr : hitsPtrVector) {
      auto const& recoDigis = hitsPtr->recoCaloDigis();
      int cryID = hitsPtr->crystalID();
      if (cryID > nCrystals * nDisks)
        continue;

      // loop over hits in the crystal
      for (auto const& rdPtr : recoDigis) {
        double VmaxHit = -1;
        int sipmID = -1;

        if (_useMeV) {
          VmaxHit = rdPtr->energyDep();
          sipmID = rdPtr->SiPMID();
        }

        else {
          auto const& rawDigiPtr = rdPtr->caloDigiPtr();
          if (rawDigiPtr.isNonnull()) {
            const std::vector<int>& waveform = rawDigiPtr->waveform();

            // evaluate baseline on first 5 samples in the waveform
            double baseline = 0;
            for (int i = 0; i < 5; ++i)
              baseline += waveform.at(i);
            baseline /= 5.0;

            // find maximum of the waveform (baseline subtracted)
            for (int i = 0; i < (int)waveform.size(); ++i) {
              if (waveform.at(i) > VmaxHit)
                VmaxHit = waveform.at(i);
            }
            if (_diagLevel > 0) {
              std::cout << "sipmID:  " << rawDigiPtr->SiPMID() << " waveform peak: " << VmaxHit
                        << " baseline: " << baseline;
            }
            VmaxHit = VmaxHit - baseline;
            sipmID = rawDigiPtr->SiPMID();
          } // caloDigi check
        }   // use MeV or ADC

        if (_diagLevel > 0)
          std::cout << " Vmax:  " << VmaxHit << std::endl;

        // cut on energy deposition
        float enecut = 9999.;
        if (cryID == 610 || cryID == 637 || cryID == 609 || cryID == 582)
          enecut = LYSOcut;
        else
          enecut = CutEnergyDep;
        if (VmaxHit > enecut && sipmID >= 0 && sipmID < nROchan) {
          float PosX = cal->mu2eToDiskFF(diskID, cal->crystal(cryID).position()).getX();
          float PosY = cal->mu2eToDiskFF(diskID, cal->crystal(cryID).position()).getY();

          PosXY.push_back({PosX, PosY});
          IDs[nhits] = sipmID;
          Vmax[nhits] = VmaxHit;
          nhits++;
        } // energy cut

      } // loop over hits
    }   // loop over crystals

    if (_diagLevel > 0)
      std::cout << "nhits: " << nhits << std::endl;

    if (nhits == 0) {
      if (_diagLevel > 0)
        std::cout << "Cluster " << iClu << " has no valid hit." << std::endl;
      continue;
    }

    // get rid of duplicate positions
    std::vector<std::pair<float, float>> fitPosXY = PosXY;
    std::sort(fitPosXY.begin(), fitPosXY.end());
    auto last = std::unique(fitPosXY.begin(), fitPosXY.end());
    fitPosXY.erase(last, fitPosXY.end());

    // find min and max x coordinate
    float min_x = fitPosXY.front().first;
    float max_x = fitPosXY.back().first;
    Dx = abs(max_x - min_x);

    if (_diagLevel > 0) {
      std::cout << "dx: " << Dx << std::endl;
    }

    // cut on ncrystalhit
    if (int(fitPosXY.size()) > CutNCryHit) {

      // prepare arrays for track fit
      std::vector<float> fitX, fitY, fitErr;
      for (auto const& p : fitPosXY) {
        fitX.push_back(p.first);
        fitY.push_back(p.second);
        fitErr.push_back(9.81);
        if (_diagLevel > 0) {
          std::cout << "x: " << p.first << " y: " << p.second << std::endl;
        }
      }

      // geometry track
      float TrkSlope = -999;
      float TrkIntercept = -999;
      float chi2norm = 999;

      if (Dx <= MaxDxVertical) { // track is vertical
        TrkSlope = 144.;
        TrkIntercept = 0.;
      } else {
        double m_init = (fitY.back() - fitY.front()) / (fitX.back() - fitX.front());
        double q_init = fitY.front() - m_init * fitX.front();

        myPoly->SetParameters(q_init, m_init);

        TGraphErrors g(fitX.size(), fitX.data(), fitY.data(), fitErr.data(), fitErr.data());
        TFitResultPtr fitresult = g.Fit(myPoly, "SQRM");

        if ((fitresult.Get() != nullptr) && (fitresult->IsValid()) && (fitresult->Ndf() > 0)) {
          TrkSlope = fitresult->Parameter(1);
          TrkIntercept = fitresult->Parameter(0);
          chi2norm = fitresult->Chi2() / fitresult->Ndf();

          if (_diagLevel > 0)
            std::cout << "Track fit done " << std::endl;
        }
      }

      // renormalize diagonal or good quality tracks
      if (Dx > MaxDxVertical) {
        for (int i = 0; i < nhits; i++) {
          if (_diagLevel > 0) {
            std::cout << "using m:" << TrkSlope << " q: " << TrkIntercept
                      << " x: " << PosXY[i].first << " y: " << PosXY[i].second
                      << " Vmax: " << Vmax[i] << std::endl;
          }
          path[i] = findpath(TrkSlope, TrkIntercept, PosXY[i].first, PosXY[i].second);
        }
      }

      // fill histograms
      for (int kk = 0; kk < nhits; kk++) {

        // vertical track: no normalization
        if (Dx <= MaxDxVertical) {
          hSiPMv[IDs[kk]]->Fill(Vmax[kk]);
          hSiPM[IDs[kk]]->Fill(Vmax[kk]);
        }

        // diagonal and general tracks
        else if ((chi2norm < CutChi2Norm) && (path[kk] > 0)) {

          // diagonal track
          if (chi2norm < 0.01)
            hSiPMd[IDs[kk]]->Fill(Vmax[kk] * cryDim / path[kk]);
          // general track
          else
            hSiPMg[IDs[kk]]->Fill(Vmax[kk] * cryDim / path[kk]);

          hSiPM[IDs[kk]]->Fill(Vmax[kk] * cryDim / path[kk]);
        }
      }
    } // end ncrystal cut
  }   // end loop iClu in ncluster
} // end analyze

void CaloCosmicEnecalib::endJob() {
  std::cout << "CaloCosmicEnecalib: entering endjob" << std::endl;

  sipmParams.assign(nROchan, LangausResult());
  sipmParamsv.assign(nROchan, LangausResult());
  sipmParamsg.assign(nROchan, LangausResult());
  sipmParamsd.assign(nROchan, LangausResult());
  int nEntries;

  for (int irou = 0; irou < nROchan; irou++) {
    nEntries = hSiPM[irou]->GetEntries();
    if (nEntries > 0) {
      if (_diagLevel > 0)
        std::cout << " >>  Fitting channel " << irou << std::endl;

      // Rebin hitsograms to help with convercence
      int rebinfactor = 1;
      if (nEntries <= lowStatThreshold)
        rebinfactor = 6;
      else if (nEntries <= midStatThreshold)
        rebinfactor = 5;
      else
        rebinfactor = 3;
      hSiPM[irou] = (TH1F*)hSiPM[irou]->Rebin(rebinfactor);

      // First attempt for langaus fit
      sipmParams[irou] = fitLangaus(hSiPM[irou]);
      sipmParamsv[irou] = fitLangaus(hSiPMv[irou]);
      sipmParamsd[irou] = fitLangaus(hSiPMd[irou]);
      sipmParamsg[irou] = fitLangaus(hSiPMg[irou]);

      hWidth->Fill(sipmParams[irou].width);
      hSigma->Fill(sipmParams[irou].sigma);
      if (sipmParamsv[irou].mpv > 0)
        MPVv->Fill(sipmParamsv[irou].mpv);
      if (sipmParamsd[irou].mpv > 0)
        MPVd->Fill(sipmParamsd[irou].mpv);
      if (sipmParamsg[irou].mpv > 0)
        MPVg->Fill(sipmParamsg[irou].mpv);
    }
  }

  // Distribution of Landau width and Gaussian sigma for all channels
  mygaus->SetParameters(hWidth->GetEntries(), hWidth->GetMean(), hWidth->GetRMS());
  hWidth->Fit("mygaus", "Q0");
  mygaus->SetParameters(hSigma->GetEntries(), hSigma->GetMean(), hSigma->GetRMS());
  hSigma->Fit("mygaus", "Q0");

  // Check fit quality
  for (int irou = 0; irou < nROchan; irou++) {
    if (hSiPM[irou]->GetEntries() > 0) {
      FitFlags flags;
      flags.fitQuality(sipmParams[irou].width, sipmParams[irou].sigma, hWidth, hSigma);

      // Refit up to 10 times if fit was not good
      int iter = 0;
      while (flags.isBad() && iter < 10) {
        if (_diagLevel > 0)
          std::cout << "Refitting channel " << irou << ", iteration " << iter << std::endl;

        sipmParams[irou] =
            fitLangaus(hSiPM[irou], 1, iter, sipmParams[irou].width, sipmParams[irou].sigma);
        flags.fitQuality(sipmParams[irou].width, sipmParams[irou].sigma, hWidth, hSigma);
        iter++;
      }
      if (sipmParams[irou].mpv > 0)
        MPV->Fill(sipmParams[irou].mpv);
    }
  }

  MPV->Fit("gaus", "Q");
  MPVv->Fit("gaus", "Q");
  MPVd->Fit("gaus", "Q");
  MPVg->Fit("gaus", "Q");

  // Write output file
  _outputfile << std::left << std::setw(6) << "#ROU" << std::setw(10) << "MPV" << std::setw(10)
              << "MPVerr" << std::setw(10) << "Width" << std::setw(10) << "Werr" << std::setw(10)
              << "Sigma" << std::setw(10) << "Serr" << std::setw(10) << "Chi2" << std::setw(8)
              << "NDF" << std::setw(8) << "Nev" << std::endl;

  for (int irou = 0; irou < nROchan; irou++) {
    const auto& res = sipmParams[irou];

    _outputfile << std::left << std::setw(6) << irou << std::fixed << std::setprecision(2)
                << std::setw(10) << res.mpv << std::setw(10) << res.mpvErr << std::setw(10)
                << res.width << std::setw(10) << res.widthErr << std::setw(10) << res.sigma
                << std::setw(10) << res.sigmaErr << std::setw(10) << res.chi2 << std::setw(8)
                << res.ndf << std::setw(8) << res.nev << std::endl;
  }

  _outputfile.close();
  delete mygaus;
  delete myPoly;
  delete mylang;

} // end endjob

double CaloCosmicEnecalib::langaus(double* x, double* par) {
  // Fit parameters:
  // par[0]=Width (scale) parameter of Landau density
  // par[1]=Most Probable (MP, location) parameter of Landau density
  // par[2]=Total area (integral -inf to inf, normalization constant)
  // par[3]=Width (sigma) of convoluted Gaussian function
  //
  // In the Landau distribution (represented by the CERNLIB approximation),
  // the maximum is located at x=-0.22278298 with the location parameter=0.
  // This shift is corrected within this function, so that the actual
  // maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  Double_t mpshift = -0.22278298;      // Landau maximum location

  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0];

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp - xlow) / np;

  // Convolution integral of Landau and Gaussian by sum
  for (i = 1.0; i <= np / 2; i++) {
    xx = xlow + (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
} // end langaus

float CaloCosmicEnecalib::findpath(float m, float q, float x, float y) {
  float xup = 0;
  float yright = 0;
  float path = 0;
  float xlow = 0;
  float yleft = 0;

  int diag = 0;
  float halfcrysize = cryDim / 2.;

  if (m != 0) {
    xup = ((y + halfcrysize) - q) / m;
    xlow = ((y - halfcrysize) - q) / m;
    yleft = m * (x - halfcrysize) + q;
    yright = m * (x + halfcrysize) + q;
  }
  if (xup < (x + halfcrysize) && xup > (x - halfcrysize)) { // hit the top of crystal
    if (xlow < (x + halfcrysize) && xlow > (x - halfcrysize)) {
      path = halfcrysize * 2 * TMath::Sqrt(1 / (m * m) + 1);
      if (diag == 1)
        std::cout << "x0,y0: (" << xup << ";" << y + halfcrysize << ") x1,y1: (" << xlow << ";"
                  << y - halfcrysize << ")" << std::endl;
    }
    if (xlow < (x - halfcrysize)) {
      path = TMath::Sqrt((xup - (x - halfcrysize)) * (xup - (x - halfcrysize)) +
                         ((y + halfcrysize) - yleft) * ((y + halfcrysize) - yleft));
      if (diag == 1) {
        std::cout << "Track exit from left" << std::endl;
        std::cout << "x0,y0: (" << xup << ";" << y + halfcrysize << ") x1,y1: (" << x - halfcrysize
                  << ";" << yleft << ")" << std::endl;
      }
    }
    if (xlow > (x + halfcrysize)) {
      path = TMath::Sqrt((xup - (x + halfcrysize)) * (xup - (x + halfcrysize)) +
                         ((y + halfcrysize) - yright) * ((y + halfcrysize) - yright));
      if (diag == 1) {
        std::cout << "x0,y0: (" << xup << ";" << y + halfcrysize << ") x1,y1: (" << x + halfcrysize
                  << ";" << yright << ")" << std::endl;
        std::cout << "Track exit from right" << std::endl;
      }
    }
  }
  // enter from left
  else if (yleft < (y + halfcrysize) && yleft > (y - halfcrysize)) {
    if (xlow < (x + halfcrysize) && xlow > (x - halfcrysize)) {
      path = TMath::Sqrt((x - halfcrysize - xlow) * (x - halfcrysize - xlow) +
                         (yleft - (y - halfcrysize)) * (yleft - (y - halfcrysize)));
    } else { // track exit from right
      path = halfcrysize * 2 * TMath::Sqrt(1 + m * m);
    }
  }

  else if (yright < (y + halfcrysize) && yright > (y - halfcrysize)) {
    if (xlow < (x + halfcrysize) && xlow > (x - halfcrysize)) {
      path = TMath::Sqrt((x + halfcrysize - xlow) * (x + halfcrysize - xlow) +
                         (yright - (y - halfcrysize)) * (yright - (y - halfcrysize)));
    }
  }

  else {
    if (diag == 1)
      std::cout << "can't find a path" << std::endl;
  }

  if (diag == 1)
    std::cout << "path is: " << path << std::endl;

  return path;
} // end findpath

CaloCosmicEnecalib::LangausResult CaloCosmicEnecalib::fitLangaus(TH1F* hist, int refit, int iter,
                                                                 float oldwidth, float oldsigma) {

  LangausResult res;

  // discarding histograms with not enough entries
  if (hist->GetEntries() < 10)
    return res;

  // set initial parameters
  int binMax = hist->GetMaximumBin();
  float maxVal = hist->GetBinContent(binMax);
  float mean = hist->GetBinCenter(binMax);
  float rms = hist->GetRMS();
  float norm = hist->GetEntries();

  // Range for Gaussian fit
  float x_lo = mean - rms / 2.;
  float x_hi = mean + rms / 2.;

  // Find FWHM
  int bin1 = binMax;
  int bin2 = binMax;
  while ((bin1 > 1) && (hist->GetBinContent(bin1) > maxVal / 2.))
    bin1--;
  while ((bin2 < hist->GetNbinsX()) && (hist->GetBinContent(bin2) > maxVal / 2.))
    bin2++;
  double FWHM = hist->GetBinCenter(bin2) - hist->GetBinCenter(bin1);

  // Perform Gaussian fit
  mygaus->SetParameters(norm, mean, FWHM / 2.);
  hist->Fit("mygaus", "Q0", "", x_lo, x_hi);
  double mpvVal = mygaus->GetParameter(1);

  // Range for Langaus fit
  double xmin = 0.;
  double xmax;
  if (!_useMeV)
    xmax = 2100.;
  else
    xmax = 100.;

  mylang->SetRange(0., xmax);

  mylang->SetParName(0, "Width");
  mylang->SetParName(1, "MPV");
  mylang->SetParName(2, "Norm");
  mylang->SetParName(3, "Sigma");

  // Initial guess for parameters in Langaus fit
  if (refit == 0)
    mylang->SetParameters(FWHM / 2., mpvVal, norm, FWHM / 2.);
  else {
    // Refit with different initial params
    float par0, par3;
    FitFlags flags;
    flags.fitQuality(oldwidth, oldsigma, hWidth, hSigma);

    // check Landau width
    if (flags.badWidthLow)
      par0 = hWidth->GetMean() + hWidth->GetRMS() * iter * 0.5;
    else if (flags.badWidthHigh)
      par0 = hWidth->GetMean() - hWidth->GetRMS() * iter * 0.5;
    else
      par0 = hWidth->GetMean();

    // check Gaus sigma
    if (flags.badSigmaLow)
      par3 = hSigma->GetMean() + hSigma->GetRMS() * iter * 0.5;
    else if (flags.badSigmaHigh)
      par3 = hSigma->GetMean() - hSigma->GetRMS() * iter * 0.5;
    else
      par3 = hSigma->GetMean();

    // Use new values as initial parameters
    mylang->SetParameters(par0, mpvVal, norm, par3);
  }

  // Limit parameters
  mylang->SetParLimits(0, xmin, hWxmax);
  mylang->SetParLimits(1, std::max(xmin, mpvVal - FWHM / 2.), std::min(xmax, mpvVal + FWHM / 2.));
  mylang->SetParLimits(3, xmin, hSxmax);

  hist->Fit("mylang", "Q", "", 0, xmax);

  if (_diagLevel > 0) {
    std::cout << "fit complete " << std::endl;
  }

  res.mpv = mylang->GetParameter(1);
  res.mpvErr = mylang->GetParError(1);
  res.width = mylang->GetParameter(0);
  res.widthErr = mylang->GetParError(0);
  res.sigma = mylang->GetParameter(3);
  res.sigmaErr = mylang->GetParError(3);
  res.chi2 = mylang->GetChisquare();
  res.ndf = mylang->GetNDF();
  res.nev = norm;

  return res;
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloCosmicEnecalib)
