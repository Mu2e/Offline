// CRV reco DQM client.
//
// Original Author: R. Mina

#include "Offline/CRVDQM/inc/CRVRecoDQM.hh"

#include "TF1.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <map>
#include <memory>

namespace {

double LandauGaussFunction(double *x, double *par)
{
    //From $ROOTSYS/tutorials/fit/langaus.C
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    constexpr Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    constexpr Double_t mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    constexpr Double_t np = 100.0;      // number of convolution steps
    constexpr Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    Double_t xx = 0.0;
    Double_t mpc = 0.0;
    Double_t fland = 0.0;
    Double_t sum = 0.0;
    Double_t xlow = 0.0, xupp = 0.0;
    Double_t step = 0.0;
    Int_t    i = 0.0;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++)
    {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);

      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}

} // end anonymous namespace for the LandauGauss function

namespace mu2e {

using namespace CRVDQMRun1;

CRVRecoDQM::CRVRecoDQM(const DQMHistSet::Config& hists) :
    DQMClient("CRVRecoDQM", kBinningVersion, hists),
    offlineToOnline_(kNOfflineChannels, -1)
{}

void CRVRecoDQM::book()
{
  DQMHistSet& h = hists();

  h_nEventsWithClusters_ = h.book1<TH1F>(
      "nEventsWithCoincidenceClusters", "Events with a coincidence cluster;;Events", kOne);
  h_coincidenceClusters_ = h.book1<TH1I>(
      "crvCoincidencesClusters", "crvCoincidenceClusters:sectorType", kSectorType);
  h_PEsVsChannel_ = h.book2<TH2F>(
      "crvPEsVsChannel",
      "PE of coincidence-cluster pulses;Offline channel (bar#times4+SiPM);PE", kOfflineChannel,
      kPE);

  // Fit results, filled once at end of job: bar charts, the MPV is the bin content.
  h_PEsMPV_ = h.bookSummary2<TH2F>(
      "crvPEsMPV", "crvPEsMPV:FEBchannel:FEB;FEB channel;FEB port", kFebChannelEdges,
      kFebPortEdges);
  h_PEsMPVROC_.resize(kNROC);
  for (int roc = 1; roc <= kNROC; ++roc) {
    h_PEsMPVROC_[roc - 1] = h.bookSummary1<TH1F>(
        Form("crvPEsMPV_ROC%d", roc),
        Form("crvPEsMPV_ROC%d;Online channel in ROC;PE MPV", roc), kRocChannelEdges);
  }
  h_PEsMPVSector_.resize(kNConfigurations);
  for (int c = 0; c < kNConfigurations; ++c) {
    for (int s = 0; s < kLayouts[c].nSectors; ++s) {
      const std::string name = "crvPEsMPV_CRVsector" + sectorTag(c, s);
      h_PEsMPVSector_[c].push_back(
          h.bookSummary1<TH1F>(name, name + ";PE MPV;Channels", kPE));
    }
  }

  // Names and titles are those of Mu2e/DQM DqmCrv_module.cc (ValCrvRecoPulse +
  // ValCrvCoincidenceCluster), so the series stays comparable across years.
  h_NPulses_ = h.book1<TH1D>("NPulses", "N Pulses", kNPulses);
  h_NPulse2_ = h.book1<TH1D>("NPulse2", "N Pulses", kNPulses2);
  h_BarIdr_ = h.book1<TH1D>("BarIdr", "RPulse Bar ID", kBarId);
  // Named on the DqmCrv "r" convention so it is not confused with the digi SiPM.
  h_SiPMr_ = h.book1<TH1D>("SiPMr", "RPulse SiPM", kSiPM);
  h_PEr_ = h.book1<TH1D>("PEr", "Fit Photoelectrons", kPulsePE);
  h_PEHeight_ = h.book1<TH1D>("PEHeight", "PE from Pulse Height", kPulsePE);
  h_PulseTime_ = h.book1<TH1D>("PulseTime", "Pulse Peak Time", kTime);
  h_PulseTime2_ = h.book1<TH1D>("PulseTime2", "Pulse Peak Time", kTime2);
  h_chi2_ = h.book1<TH1D>("chi2", "Pulse fit chi2", kChi2);
  h_logchi2_ = h.book1<TH1D>("logchi2", "log10 Pulse fit chi2", kLogChi2);
  h_LeadingTime_ = h.book1<TH1D>("LeadingTime", "Leading Edge Time", kTime);
  h_LeadingTime2_ = h.book1<TH1D>("LeadingTime2", "Leading Edge Time", kTime2);

  h_NClus_ = h.book1<TH1D>("NClus", "N Clusters", kNClusters);
  h_NPc_ = h.book1<TH1D>("NPc", "N Pulse", kNClusters);
  h_PEc_ = h.book1<TH1D>("PEc", "clus PE", kClusterPE);
  h_tc_ = h.book1<TH1D>("tc", "clus start time", kTime);
  h_t2c_ = h.book1<TH1D>("t2c", "clus start time", kTime2);

  static const char* kCoordinate[3] = {"X", "Y", "Z"};
  h_position_.resize(kNConfigurations);
  for (int c = 0; c < kNConfigurations; ++c) {
    for (int i = 0; i < 3; ++i) {
      h_position_[c].push_back(h.book1<TH1D>(
          Form("%s_%s", kCoordinate[i], kLayouts[c].name),
          Form("clus %s (%s);%s [mm];Clusters", kCoordinate[i], kLayouts[c].name,
               kCoordinate[i]),
          positionAxis(c, i)));
    }
  }
}

TH1F* CRVRecoDQM::PEsMPVROC(int roc) const
{
  if (roc < 1 || static_cast<std::size_t>(roc) > h_PEsMPVROC_.size()) {
    return nullptr;
  }
  return h_PEsMPVROC_[roc - 1];
}

TH1F* CRVRecoDQM::PEsMPVSector(int configuration, int sector) const
{
  if (configuration < 0 || static_cast<std::size_t>(configuration) >= h_PEsMPVSector_.size() ||
      sector < 0 || static_cast<std::size_t>(sector) >= h_PEsMPVSector_[configuration].size()) {
    return nullptr;
  }
  return h_PEsMPVSector_[configuration][sector];
}

TH1D* CRVRecoDQM::position(int configuration, int coordinate) const
{
  if (configuration < 0 || static_cast<std::size_t>(configuration) >= h_position_.size()) {
    return nullptr;
  }
  return h_position_[configuration][coordinate];
}

void CRVRecoDQM::SetConfiguration(int configuration, const std::vector<int>& channelToSector)
{
  configuration_ = configuration;
  channelToSector_ = channelToSector;
}

void CRVRecoDQM::Fill(const CrvCoincidenceClusterCollection& clusters)
{
  fillClusters(clusters);
}

void CRVRecoDQM::Fill(const CrvCoincidenceClusterCollection& clusters,
                      const CrvRecoPulseCollection& recoPulses)
{
  fillClusters(clusters);
  fillRecoPulses(recoPulses);
}

void CRVRecoDQM::fillClusters(const CrvCoincidenceClusterCollection& clusters)
{
  if (!booked()) {
    return;
  }
  // No CrvStatus on this tier, so the window clock counts events.
  beginEvent();

  if (!clusters.empty()) {
    ++nEventsWithClusters_;
    h_nEventsWithClusters_.Fill(1.);
  }
  //an entry every event, zero included: that is what makes NClus a rate
  h_NClus_.Fill(clusters.size());

  if (configuration_ < 0 && !clusters.empty()) {
    diag().Count("noConfiguration",
                 "no CRV configuration was set; cluster positions are not filled.");
  }

  for (const auto& cluster : clusters) {
    ++nClusters_;
    h_coincidenceClusters_.Fill(cluster.GetCrvSectorType());
    h_NPc_.Fill(cluster.GetCrvRecoPulses().size());
    h_PEc_.Fill(cluster.GetPEs());
    h_tc_.Fill(cluster.GetStartTime());
    h_t2c_.Fill(cluster.GetStartTime());
    if (configuration_ >= 0) {
      const auto& pos = cluster.GetAvgHitPos();
      bool outside = false;
      for (int i = 0; i < 3; ++i) {
        const DQMAxis a = positionAxis(configuration_, i);
        outside = outside || pos[i] < a.lo || pos[i] >= a.hi;
        h_position_[configuration_][i].Fill(pos[i]);
      }
      if (outside) {
        diag().Count("clusterOutsideEnvelope",
                     "a cluster position is outside its configuration's X/Y/Z axes "
                     "(it is in under/overflow).");
      }
    }

    for (const auto& pulsePtr : cluster.GetCrvRecoPulses()) {
      if (!pulsePtr) {
        diag().Count("nullRecoPulsePtr",
                     "a coincidence cluster holds a null CrvRecoPulse Ptr, usually because "
                     "the CrvRecoPulse product is not in the event; the pulse is left out "
                     "of crvPEsVsChannel.");
        continue;
      }
      const CrvRecoPulse& pulse = *pulsePtr;
      ++nRecoPulses_;

      const int channel = pulse.GetScintillatorBarIndex().asInt() *
                              static_cast<int>(CRVId::nChanPerBar) +
                          pulse.GetSiPMNumber();
      if (channel < 0 || channel >= kNOfflineChannels) {
        diag().Count("offlineChannelOutOfRange",
                     "a reco pulse's offline channel is outside CRVId; it is left out of "
                     "crvPEsVsChannel.");
        continue;
      }
      h_PEsVsChannel_.Fill(channel, pulse.GetPEs());

      // ROC/FEB/FEBchannel come off the digi the pulse was made from.
      const int roc = static_cast<int>(pulse.GetROC());
      const int feb = static_cast<int>(pulse.GetFEB());
      const int febChannel = static_cast<int>(pulse.GetFEBchannel());
      if (onlineIdInRange(roc, feb, febChannel)) {
        offlineToOnline_[channel] = onlineChannel(roc, feb, febChannel);
      } else {
        diag().Count("onlineIdOutOfRange",
                     "a reco pulse's ROC/FEB/channel is outside CRVId; the online MPV "
                     "maps leave its channel out.");
      }
    }
  }
}

void CRVRecoDQM::fillRecoPulses(const CrvRecoPulseCollection& recoPulses)
{
  if (!booked()) {
    return;
  }
  h_NPulses_.Fill(recoPulses.size());
  h_NPulse2_.Fill(recoPulses.size());

  for (const auto& pulse : recoPulses) {
    h_BarIdr_.Fill(pulse.GetScintillatorBarIndex().asInt());
    h_SiPMr_.Fill(pulse.GetSiPMNumber());
    h_PEr_.Fill(pulse.GetPEs());
    h_PEHeight_.Fill(pulse.GetPEsPulseHeight());
    h_PulseTime_.Fill(pulse.GetPulseTime());
    h_PulseTime2_.Fill(pulse.GetPulseTime());
    h_chi2_.Fill(pulse.GetPulseFitChi2());
    //a chi2 of zero means the fit did not run, and log10 of it is not a number
    if (pulse.GetPulseFitChi2() > 0.0) {
      h_logchi2_.Fill(std::log10(pulse.GetPulseFitChi2()));
    }
    h_LeadingTime_.Fill(pulse.GetLEtime());
    h_LeadingTime2_.Fill(pulse.GetLEtime());
  }
}

CRVRecoDQM::LandauGaussResult CRVRecoDQM::FitSpectrum(TH1* h)
{
  LandauGaussResult result;
  if (h == nullptr) {
    return result;
  }

  std::multimap<float, float> bins;  //binContent,binCenter
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    if (h->GetBinCenter(i) < kFitMinBinCenter) continue;
    bins.emplace(h->GetBinContent(i), h->GetBinCenter(i));  //ordered from smallest to largest bin entries
  }
  if (bins.size() < 4) return result;
  if (bins.rbegin()->first < kFitMinPeakEntries) return result;  //low statistics

  int nBins = 0;
  float binSum = 0;
  for (auto bin = bins.rbegin(); bin != bins.rend(); ++bin) {
    nBins++;
    binSum += bin->second;
    if (nBins == 4) break;
  }
  float maxX = binSum / 4;
  float fitRangeStart = kFitRangeStart * maxX;
  float fitRangeEnd = kFitRangeEnd * maxX;
  if (maxX < kFitMinPE) maxX = kFitMinPE;
  if (fitRangeStart < kFitMinPE) fitRangeStart = kFitMinPE;

  //Parameters
  Double_t startValues[4], parLimitsLow[4], parLimitsHigh[4];
  //Most probable value
  startValues[1] = maxX;
  parLimitsLow[1] = fitRangeStart;
  parLimitsHigh[1] = fitRangeEnd;
  //Area
  startValues[2] = h->Integral(h->FindBin(fitRangeStart), h->FindBin(fitRangeEnd));
  parLimitsLow[2] = 0.01 * startValues[2];
  parLimitsHigh[2] = 100 * startValues[2];
  //Other parameters
  startValues[0] = 5.0;    startValues[3] = 10.0;
  parLimitsLow[0] = 2.0;   parLimitsLow[3] = 2.0;
  parLimitsHigh[0] = 15.0; parLimitsHigh[3] = 20.0;

  TF1 fit("LandauGauss", LandauGaussFunction, fitRangeStart, fitRangeEnd, 4);
  fit.SetParameters(startValues);
  fit.SetParNames("Width", "MP", "Area", "GSigma");
  for (int i = 0; i < 4; i++) fit.SetParLimits(i, parLimitsLow[i], parLimitsHigh[i]);
  TFitResultPtr fr = h->Fit(&fit, "LQRSN");

  const float mpv = fit.GetMaximumX();
  //a fit that ran but landed on the range edge did not find a peak
  if (mpv == fitRangeStart) return result;

  result.mpv = mpv;
  result.fitted = true;
  if (fr.Get() != nullptr && fr->Ndf() > 0) {
    result.chi2PerNdf = fr->Chi2() / fr->Ndf();
  }
  return result;
}

std::vector<CRVRecoDQM::LandauGaussResult> CRVRecoDQM::FitPEsVsChannel(const TH2& peVsChannel)
{
  const int nChannels = peVsChannel.GetNbinsX();
  std::vector<LandauGaussResult> results(nChannels);
  const TAxis* y = peVsChannel.GetYaxis();
  auto spectrum = std::make_unique<TH1D>("crvPEsSpectrum", "", y->GetNbins(), y->GetXmin(),
                                         y->GetXmax());
  spectrum->SetDirectory(nullptr);
  for (int ix = 1; ix <= nChannels; ++ix) {
    double total = 0.;
    for (int iy = 1; iy <= y->GetNbins(); ++iy) {
      const double c = peVsChannel.GetBinContent(ix, iy);
      spectrum->SetBinContent(iy, c);
      total += c;
    }
    if (total > 0.) {
      results[ix - 1] = FitSpectrum(spectrum.get());
    }
  }
  return results;
}

void CRVRecoDQM::endJob()
{
  const std::vector<LandauGaussResult> results = FitPEsVsChannel(*h_PEsVsChannel_.get());

  for (std::size_t channel = 0; channel < results.size(); ++channel) {
    const int online = offlineToOnline_[channel];
    if (online < 0) {
      continue;
    }
    ++nFits_;
    if (results[channel].fitted) {
      ++nFitsSucceeded_;
    }
    const int port = online / kNChanPerFEB;
    const int febChannel = online % kNChanPerFEB;
    const int roc = port / kNFebPerROC + 1;
    const int feb = port % kNFebPerROC + 1;
    h_PEsMPVROC_[roc - 1].Fill(rocChannel(feb, febChannel), results[channel].mpv);
    h_PEsMPV_.Fill(febChannel, port, results[channel].mpv);
  }

  if (configuration_ < 0) {
    return;
  }
  //an unfittable channel enters as zero, which is what makes it visible
  const auto& family = h_PEsMPVSector_[configuration_];
  for (std::size_t channel = 0; channel < channelToSector_.size() && channel < results.size();
       ++channel) {
    const int sector = channelToSector_[channel];
    if (sector >= 0 && static_cast<std::size_t>(sector) < family.size()) {
      family[sector].Fill(results[channel].mpv);
    }
  }
}

} // namespace mu2e
