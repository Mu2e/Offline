// Standalone CRV reco DQM helper.
//
// Original Author: R. Mina

#include "Offline/DQMHelpers/inc/CRVRecoDQM.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TF1.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TString.h"

#include <cmath>
#include <map>

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

// Bins below this PE value are dark noise and the reco threshold turn-on, and
// would drag the peak estimate down.
constexpr double kFitMinBinCenter = 10.0;
// Fewer than this in the tallest bin is not a spectrum worth fitting.
constexpr double kFitMinPeakEntries = 20.0;

} // end anonymous namespace for the LandauGauss function

namespace mu2e {

bool CRVRecoDQM::onlineIdInRange(int roc, int feb, int febChannel)
{
  return roc >= 1 && roc <= kNROC && feb >= 1 && feb <= kNFebPerROC &&
         febChannel >= 0 && febChannel < kNChanPerFEB;
}

int CRVRecoDQM::rocChannel(int feb, int febChannel)
{
  return (feb - 1) * kNChanPerFEB + febChannel;
}

int CRVRecoDQM::febPort(int roc, int feb)
{
  return (roc - 1) * kNFebPerROC + feb - 1;
}

int CRVRecoDQM::onlineChannelIndex(int roc, int feb, int febChannel)
{
  return (roc - 1) * kNChannelsPerROC + rocChannel(feb, febChannel);
}

namespace {
// ROOT quietly makes an unusable axis from a non-positive bin count or an
// inverted range, and the plot then looks broken rather than misconfigured.
void checkAxis(const char* what, int& nBins, double lo, double hi)
{
  if (nBins < 1) {
    mf::LogWarning("CRVRecoDQM")
        << what << " bin count " << nBins << " is not positive; using 1.";
    nBins = 1;
  }
  if (!(hi > lo)) {
    mf::LogWarning("CRVRecoDQM")
        << what << " range [" << lo << ", " << hi << "] is empty or inverted. "
        << "ROOT will not make a usable axis from it; the plot will be unreadable.";
  }
}
} // namespace

CRVRecoDQM::CRVRecoDQM(const Config& config) : config_(config)
{
  checkAxis("PE", config_.nBinsPEs, config_.minPEs, config_.maxPEs);
  if (config_.nSectorTypeBins < 1) {
    config_.nSectorTypeBins = 1;
  }
  if (config_.fillInclusive) {
    checkAxis("time", config_.nBinsTime, config_.minTime, config_.maxTime);
    checkAxis("long time", config_.nBinsTime2, config_.minTime2, config_.maxTime2);
    checkAxis("X", config_.nBinsPos, config_.minX, config_.maxX);
    checkAxis("Y", config_.nBinsPos, config_.minY, config_.maxY);
    checkAxis("Z", config_.nBinsPos, config_.minZ, config_.maxZ);
  }
}

void CRVRecoDQM::Book(art::TFileDirectory dir)
{
  dir_ = dir;

  h_nEvents_ = dir.make<TH1F>("nEvents", "Events processed;;Events", 1, 0.5, 1.5);
  h_nEventsWithClusters_ =
      dir.make<TH1F>("nEventsWithCoincidenceClusters",
                     "Events with a coincidence cluster;;Events", 1, 0.5, 1.5);

  h_coincidenceClusters_ =
      dir.make<TH1I>("crvCoincidencesClusters", "crvCoincidenceClusters:sectorType",
                     config_.nSectorTypeBins, 0, config_.nSectorTypeBins);

  h_PEsMPVROC_.assign(kNROC, nullptr);
  for (int roc = 1; roc <= kNROC; ++roc) {
    h_PEsMPVROC_[roc - 1] = dir.make<TH1F>(
        Form("crvPEsMPV_ROC%d", roc),
        Form("crvPEsMPV_ROC%d;Online channel in ROC;PE MPV", roc),
        kNChannelsPerROC, 0, kNChannelsPerROC);
  }

  h_PEsMPV_ = dir.make<TH2F>("crvPEsMPV", "crvPEsMPV:FEBchannel:FEB;FEB channel;FEB port",
                             kNChanPerFEB, 0, kNChanPerFEB,
                             kNFebPorts, 0, kNFebPorts);

  // The spectra are the fit input, not a DQM product, and there are ~50k of
  // them across all of CRVId. They are booked on first fill.
  h_PEsOffline_.assign(CRVId::nChannels, nullptr);
  h_PEsOnline_.assign(kNOnlineChannels, nullptr);
  if (config_.writePerChannelPE) {
    spectraDir_ = dir.mkdir("perChannelPE");
  }

  if (config_.fillInclusive) {
    bookInclusive(dir);
  }

  booked_ = true;
}

// Names and titles are those of Mu2e/DQM DqmCrv_module.cc, which are in turn
// ValCrvRecoPulse.cc + ValCrvCoincidenceCluster.cc. Keep them: the point of
// these plots is a series that stays comparable across years. Only the axes
// that depend on the detector's position and readout window are configurable;
// see Config.
void CRVRecoDQM::bookInclusive(art::TFileDirectory& dir)
{
  h_NPulses_ = dir.make<TH1D>("NPulses", "N Pulses", 101, -0.5, 100.5);
  h_NPulse2_ = dir.make<TH1D>("NPulse2", "N Pulses", 100, -0.5, 2000.0);
  h_BarIdr_ = dir.make<TH1D>("BarIdr", "RPulse Bar ID", 200, -0.5, 5503.5);
  // ValCrvRecoPulse has this and the DqmCrv copy dropped it. The digi SiPM
  // hist tests the readout; this one tests reconstruction per readout end, and
  // the ratio isolates an end where waveforms arrive but do not fit. Named
  // SiPMr on the DqmCrv "r" convention so it cannot be confused with the digi
  // SiPM hist that CRVDigiDQM books.
  h_SiPMr_ = dir.make<TH1D>("SiPMr", "RPulse SiPM", 4, -0.5, 3.5);
  h_PEr_ = dir.make<TH1D>("PEr", "Fit Photoelectrons", 100, 0.0, 400.0);
  h_PEHeight_ =
      dir.make<TH1D>("PEHeight", "PE from Pulse Height", 100, 0.0, 400.0);
  h_PulseTime_ = dir.make<TH1D>("PulseTime", "Pulse Peak Time",
                                config_.nBinsTime, config_.minTime, config_.maxTime);
  h_PulseTime2_ = dir.make<TH1D>("PulseTime2", "Pulse Peak Time",
                                 config_.nBinsTime2, config_.minTime2, config_.maxTime2);
  h_chi2_ = dir.make<TH1D>("chi2", "Pulse fit chi2", 100, 0.0, 20.0);
  h_logchi2_ = dir.make<TH1D>("logchi2", "log10 Pulse fit chi2", 100, -3.0, 5.0);
  h_LeadingTime_ = dir.make<TH1D>("LeadingTime", "Leading Edge Time",
                                  config_.nBinsTime, config_.minTime, config_.maxTime);
  h_LeadingTime2_ = dir.make<TH1D>("LeadingTime2", "Leading Edge Time",
                                   config_.nBinsTime2, config_.minTime2, config_.maxTime2);

  h_NClus_ = dir.make<TH1D>("NClus", "N Clusters", 101, -0.5, 100.5);
  h_NPc_ = dir.make<TH1D>("NPc", "N Pulse", 101, -0.5, 100.5);
  h_PEc_ = dir.make<TH1D>("PEc", "clus PE", 100, 0.0, 2000.0);
  h_tc_ = dir.make<TH1D>("tc", "clus start time",
                         config_.nBinsTime, config_.minTime, config_.maxTime);
  h_t2c_ = dir.make<TH1D>("t2c", "clus start time",
                          config_.nBinsTime2, config_.minTime2, config_.maxTime2);
  // X/Y/Z are not booked here: see BookPositionAxes.

  // SecType is deliberately absent: crvCoincidencesClusters already is it, and
  // the DQM audit asks the merge to end with one of the two, not both.
}

TH1F* CRVRecoDQM::spectrum(std::vector<TH1F*>& hists, std::size_t index,
                           const char* namePrefix)
{
  if (index >= hists.size()) {
    return nullptr;
  }
  if (hists[index] != nullptr) {
    return hists[index];
  }

  const char* name = Form("%s%zu", namePrefix, index);
  if (spectraDir_) {
    hists[index] = spectraDir_->make<TH1F>(name, name, config_.nBinsPEs,
                                           config_.minPEs, config_.maxPEs);
    return hists[index];
  }
  // Belongs to no TDirectory, so nothing writes or deletes it behind our back.
  auto h = std::make_unique<TH1F>(name, name, config_.nBinsPEs, config_.minPEs,
                                  config_.maxPEs);
  h->SetDirectory(nullptr);
  hists[index] = h.get();
  ownedSpectra_.push_back(std::move(h));
  return hists[index];
}

void CRVRecoDQM::BookPositionAxes(const std::vector<double>& minPoint,
                                  const std::vector<double>& maxPoint,
                                  double margin)
{
  if (minPoint.size() < 3 || maxPoint.size() < 3) {
    mf::LogWarning("CRVRecoDQM")
        << "BookPositionAxes got a " << minPoint.size() << "/" << maxPoint.size()
        << "-component envelope instead of 3; using the configured X/Y/Z ranges "
        << "instead, which may not suit this geometry.";
    BookPositionAxes();
    return;
  }
  // A geometry with no CRV sectors yields a zero-extent box rather than an
  // error, which would make three degenerate axes.
  if (!(maxPoint[0] > minPoint[0] && maxPoint[1] > minPoint[1] &&
        maxPoint[2] > minPoint[2])) {
    mf::LogWarning("CRVRecoDQM")
        << "BookPositionAxes got an empty CRV envelope (x " << minPoint[0]
        << ".." << maxPoint[0] << ", y " << minPoint[1] << ".." << maxPoint[1]
        << ", z " << minPoint[2] << ".." << maxPoint[2]
        << "); using the configured X/Y/Z ranges instead.";
    BookPositionAxes();
    return;
  }
  // The envelope comes from the counters; a cluster position derived from the
  // two-sided time difference carries a resolution that can put it just outside.
  double lo[3], hi[3];
  for (int i = 0; i < 3; ++i) {
    const double pad = std::abs(maxPoint[i] - minPoint[i]) * margin;
    lo[i] = minPoint[i] - pad;
    hi[i] = maxPoint[i] + pad;
  }
  bookPositionHists(lo, hi);
}

void CRVRecoDQM::BookPositionAxes()
{
  const double lo[3] = {config_.minX, config_.minY, config_.minZ};
  const double hi[3] = {config_.maxX, config_.maxY, config_.maxZ};
  bookPositionHists(lo, hi);
}

void CRVRecoDQM::bookPositionHists(const double lo[3], const double hi[3])
{
  //first caller wins, so an injected envelope beats the lazy Config fallback
  if (!dir_ || !config_.fillInclusive || h_X_ != nullptr) {
    return;
  }
  checkAxis("X", config_.nBinsPos, lo[0], hi[0]);
  checkAxis("Y", config_.nBinsPos, lo[1], hi[1]);
  checkAxis("Z", config_.nBinsPos, lo[2], hi[2]);
  h_X_ = dir_->make<TH1D>("X", "clus X", config_.nBinsPos, lo[0], hi[0]);
  h_Y_ = dir_->make<TH1D>("Y", "clus Y", config_.nBinsPos, lo[1], hi[1]);
  h_Z_ = dir_->make<TH1D>("Z", "clus Z", config_.nBinsPos, lo[2], hi[2]);
}

void CRVRecoDQM::ensurePositionAxes()
{
  if (h_X_ == nullptr) {
    BookPositionAxes();
  }
}

void CRVRecoDQM::BookSectorMPV(const std::vector<std::string>& sectorNames,
                               const std::vector<int>& channelToSector)
{
  if (!dir_ || !h_PEsMPVSector_.empty()) {
    return;
  }
  channelToSector_ = channelToSector;
  h_PEsMPVSector_.reserve(sectorNames.size());
  for (const auto& sector : sectorNames) {
    const std::string name = "crvPEsMPV_CRVsector" + sector;
    h_PEsMPVSector_.push_back(dir_->make<TH1F>(name.c_str(), name.c_str(),
                                               config_.nBinsPEs, config_.minPEs,
                                               config_.maxPEs));
  }
}

TH1F* CRVRecoDQM::PEsOffline(std::size_t channel) const
{
  return channel < h_PEsOffline_.size() ? h_PEsOffline_[channel] : nullptr;
}

TH1F* CRVRecoDQM::PEsOnline(std::size_t onlineChannel) const
{
  return onlineChannel < h_PEsOnline_.size() ? h_PEsOnline_[onlineChannel]
                                             : nullptr;
}

void CRVRecoDQM::Fill(const CrvCoincidenceClusterCollection& clusters)
{
  fillClusters(clusters);
}

void CRVRecoDQM::Fill(const CrvCoincidenceClusterCollection& clusters,
                      const CrvRecoPulseCollection& recoPulses)
{
  fillClusters(clusters);
  fillInclusiveRecoPulses(recoPulses);
}

void CRVRecoDQM::fillClusters(const CrvCoincidenceClusterCollection& clusters)
{
  if (!booked_) {
    return;
  }

  ++nEvents_;
  h_nEvents_->Fill(1.f);
  if (!clusters.empty()) {
    ++nEventsWithClusters_;
    h_nEventsWithClusters_->Fill(1.f);
  }
  //an entry every event, zero included: that is what makes NClus a rate
  if (config_.fillInclusive) {
    ensurePositionAxes();  //no-op once booked, from geometry or from Config
    h_NClus_->Fill(clusters.size());
  }

  for (const auto& cluster : clusters) {
    ++nClusters_;
    h_coincidenceClusters_->Fill(cluster.GetCrvSectorType());

    if (config_.fillInclusive) {
      h_NPc_->Fill(cluster.GetCrvRecoPulses().size());
      h_PEc_->Fill(cluster.GetPEs());
      h_tc_->Fill(cluster.GetStartTime());
      h_t2c_->Fill(cluster.GetStartTime());
      if (h_X_) {
        h_X_->Fill(cluster.GetAvgHitPos().x());
        h_Y_->Fill(cluster.GetAvgHitPos().y());
        h_Z_->Fill(cluster.GetAvgHitPos().z());
      }
    }

    for (const auto& pulsePtr : cluster.GetCrvRecoPulses()) {
      if (!pulsePtr) {
        ++nNullRecoPulsePtrs_;
        if (!warnedNullRecoPulsePtr_) {
          warnedNullRecoPulsePtr_ = true;
          mf::LogWarning("CRVRecoDQM")
              << "coincidence cluster holds a null CrvRecoPulse Ptr; the pulse "
              << "is omitted from the PE spectra. This usually means the "
              << "CrvRecoPulse product the clusters point at is not in the "
              << "event. Reported once per job.";
        }
        continue;
      }

      const CrvRecoPulse& pulse = *pulsePtr;
      ++nRecoPulses_;
      const float PEs = pulse.GetPEs();

      const std::size_t channel =
          static_cast<std::size_t>(pulse.GetScintillatorBarIndex().asUint()) *
              CRVId::nChanPerBar +
          static_cast<std::size_t>(pulse.GetSiPMNumber());
      if (channel < h_PEsOffline_.size()) {
        spectrum(h_PEsOffline_, channel, "crvPEs_channel")->Fill(PEs);
      } else {
        ++nOfflineChannelOutOfRange_;
        if (!warnedOfflineChannelOutOfRange_) {
          warnedOfflineChannelOutOfRange_ = true;
          mf::LogWarning("CRVRecoDQM")
              << "reco pulse on offline channel " << channel
              << " is past the " << h_PEsOffline_.size()
              << "-channel CRVId range; it is omitted from crvPEsMPV_CRVsector*. "
              << "Reported once per job.";
        }
      }

      // ROC/FEB/FEBchannel come off the digi the pulse was made from, so the
      // online view needs no CRVOrdinal channel map here.
      const int roc = static_cast<int>(pulse.GetROC());
      const int feb = static_cast<int>(pulse.GetFEB());
      const int febChannel = static_cast<int>(pulse.GetFEBchannel());
      if (onlineIdInRange(roc, feb, febChannel)) {
        const std::size_t onlineChannel = onlineChannelIndex(roc, feb, febChannel);
        spectrum(h_PEsOnline_, onlineChannel, "crvPEsROC_channel")->Fill(PEs);
      } else {
        ++nOnlineIdOutOfRange_;
        if (!warnedOnlineIdOutOfRange_) {
          warnedOnlineIdOutOfRange_ = true;
          mf::LogWarning("CRVRecoDQM")
              << "reco pulse from ROC " << roc << " FEB " << feb << " channel "
              << febChannel << " is outside CRVId range (ROC 1-" << kNROC
              << ", FEB 1-" << kNFebPerROC << ", channel 0-"
              << (kNChanPerFEB - 1) << "). It is omitted from crvPEsMPV_ROC* "
              << "and the 2D crvPEsMPV; the per-sector MPV still gets it. "
              << "nOnlineIdOutOfRange() counts them. Reported once per job.";
        }
      }
    }
  }
}

void CRVRecoDQM::fillInclusiveRecoPulses(const CrvRecoPulseCollection& recoPulses)
{
  if (!booked_ || !config_.fillInclusive) {
    return;
  }

  nInclusiveRecoPulses_ += recoPulses.size();
  h_NPulses_->Fill(recoPulses.size());
  h_NPulse2_->Fill(recoPulses.size());

  for (const auto& pulse : recoPulses) {
    h_BarIdr_->Fill(pulse.GetScintillatorBarIndex().asInt());
    h_SiPMr_->Fill(pulse.GetSiPMNumber());
    h_PEr_->Fill(pulse.GetPEs());
    h_PEHeight_->Fill(pulse.GetPEsPulseHeight());
    h_PulseTime_->Fill(pulse.GetPulseTime());
    h_PulseTime2_->Fill(pulse.GetPulseTime());
    h_chi2_->Fill(pulse.GetPulseFitChi2());
    //a chi2 of zero means the fit did not run, and log10 of it is not a number
    if (pulse.GetPulseFitChi2() > 0.0) {
      h_logchi2_->Fill(std::log10(pulse.GetPulseFitChi2()));
    }
    h_LeadingTime_->Fill(pulse.GetLEtime());
    h_LeadingTime2_->Fill(pulse.GetLEtime());
  }
}

CRVRecoDQM::LandauGaussResult CRVRecoDQM::fitChannel(TH1F* h)
{
  LandauGaussResult result;
  if (h == nullptr) {
    return result;
  }
  ++nFits_;

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
  float fitRangeStart = config_.PEfitRangeStart * maxX;  //0.6 @ 24
  float fitRangeEnd = config_.PEfitRangeEnd * maxX;
  if (maxX < config_.PEstart) maxX = config_.PEstart;
  if (fitRangeStart < config_.PEstart) fitRangeStart = config_.PEstart;

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
  parLimitsHigh[0] = 15.0; parLimitsHigh[3] = 20.0;  //7 and 15 @ 21  //6 and 13 @ 23

  TF1 fit("LandauGauss", LandauGaussFunction, fitRangeStart, fitRangeEnd, 4);
  fit.SetParameters(startValues);
  fit.SetLineColor(kRed);
  fit.SetParNames("Width", "MP", "Area", "GSigma");
  for (int i = 0; i < 4; i++) fit.SetParLimits(i, parLimitsLow[i], parLimitsHigh[i]);
  TFitResultPtr fr = h->Fit(&fit, "LQRS");

  const float mpv = fit.GetMaximumX();
  //a fit that ran but landed on the range edge did not find a peak
  if (mpv == fitRangeStart) return result;

  result.mpv = mpv;
  result.fitted = true;
  ++nFitsSucceeded_;
  if (fr.Get() != nullptr && fr->Ndf() > 0) {
    result.chi2PerNdf = fr->Chi2() / fr->Ndf();
    fitChi2Sum_ += result.chi2PerNdf;
    ++fitChi2N_;
  }
  return result;
}

void CRVRecoDQM::fitSectorMPV()
{
  if (h_PEsMPVSector_.empty()) {
    return;
  }
  for (std::size_t channel = 0; channel < channelToSector_.size(); ++channel) {
    const int sector = channelToSector_[channel];
    if (sector < 0 ||
        static_cast<std::size_t>(sector) >= h_PEsMPVSector_.size()) {
      continue;
    }
    const LandauGaussResult result = fitChannel(PEsOffline(channel));
    //an unfittable channel enters as zero, which is what makes it visible
    h_PEsMPVSector_[sector]->Fill(result.mpv);
  }
}

void CRVRecoDQM::fitOnlineMPV()
{
  for (int roc = 1; roc <= kNROC; ++roc) {
    for (int feb = 1; feb <= kNFebPerROC; ++feb) {
      for (int febChannel = 0; febChannel < kNChanPerFEB; ++febChannel) {
        const LandauGaussResult result =
            fitChannel(PEsOnline(onlineChannelIndex(roc, feb, febChannel)));
        //bar charts, not distributions: the MPV is the bin content
        h_PEsMPVROC_[roc - 1]->Fill(rocChannel(feb, febChannel), result.mpv);
        h_PEsMPV_->Fill(febChannel, febPort(roc, feb), result.mpv);
      }
    }
  }
}

double CRVRecoDQM::meanFitChi2PerNdf() const
{
  if (fitChi2N_ == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return fitChi2Sum_ / static_cast<double>(fitChi2N_);
}

void CRVRecoDQM::WriteGraphs()
{
  if (!booked_) {
    return;
  }
  fitSectorMPV();
  fitOnlineMPV();
}

} // namespace mu2e
