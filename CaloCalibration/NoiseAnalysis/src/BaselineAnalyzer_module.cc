#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "artdaq-core/Data/Fragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>

#include "Offline/CaloVisualizer/inc/THMu2eCaloDisk.hh"
#include "Offline/DAQ/inc/CaloDAQUtilities.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "artdaq-core-mu2e/Data/EventHeader.hh"
#include "artdaq-core-mu2e/Overlays/DTCEventFragment.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/CalorimeterDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include "cetlib_except/exception.h"

//-- insert calls to proditions ..for calodmap-----
#include "Offline/CaloConditions/inc/CalDAQMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
//-------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TEllipse.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "art_root_io/TFileService.h"
#include <sys/stat.h>

namespace mu2e {
class BaselineAnalyzer : public art::EDAnalyzer {
public:
  // clang-format off
    struct Config {
      fhicl::Atom<std::string> caloDigiTag {fhicl::Name("caloDigiTag" ) , fhicl::Comment("caloDigiTag"), ""};
      fhicl::Atom<int> verbosity {fhicl::Name("verbosity" ) , fhicl::Comment("Verbosity [0-2]"), 1};
      fhicl::Atom<bool> writeCSV {fhicl::Name("writeCSV" ) , fhicl::Comment("Write CSV file with thresholds"), false};
      fhicl::Atom<std::string> CSVfilename {fhicl::Name("CSVfilename" ) , fhicl::Comment("CSV file to write thresholds into"), ""};
      fhicl::Atom<bool> writePDF {fhicl::Name("writePDF" ) , fhicl::Comment("Write PDF report"), false};
      fhicl::Atom<std::string> PDFfilename {fhicl::Name("PDFfilename" ) , fhicl::Comment("PDF report file name"), ""};
      fhicl::Sequence<int> thresholdOffset {fhicl::Name("thresholdOffset" ) , fhicl::Comment("Offsets with respect to gaussian mean (one set of outputs per combination of offset, offsetPin and Nsigma)"), std::vector<int>{100}};
      fhicl::Sequence<int> thresholdOffsetPin {fhicl::Name("thresholdOffsetPin" ) , fhicl::Comment("Offsets with respect to gaussian mean (pin diodes)"), std::vector<int>{100}};
      fhicl::Sequence<double> Nsigma {fhicl::Name("Nsigma" ) , fhicl::Comment("Threshold = mean + Nsigma*RMS + offset (0 -> mean + offset)"), std::vector<double>{0.}};
      fhicl::Atom<int> hotStdDev {fhicl::Name("hotStdDev" ) , fhicl::Comment("StdDev limit for hot channels"), 6};
      fhicl::Atom<int> coldStdDev {fhicl::Name("coldStdDev" ) , fhicl::Comment("StdDev limit for cold channels"), 3};
      fhicl::Atom<int> maxSamples {fhicl::Name("maxSamples" ) , fhicl::Comment("Maximum number of samples to use for baseline calculation. -1 for all samples"), -1};
      fhicl::Atom<int> minSubrunSamples {fhicl::Name("minSubrunSamples" ) , fhicl::Comment("Minimum number of samples of a channel in a subrun to compute its average baseline"), 100};
      fhicl::Atom<int> maxPeakToPeak {fhicl::Name("maxPeakToPeak" ) , fhicl::Comment("Waveforms with max-min >= this value [ADC] are discarded (e.g. cosmic ray signals)"), 50};
      fhicl::Atom<double> minimumStdDev {fhicl::Name("minimumStdDev" ) , fhicl::Comment("Minimum StdDev"), 3.};
      fhicl::Atom<int> minimumThreshold {fhicl::Name("minimumThreshold" ) , fhicl::Comment("Minimum threshold above the baseline [ADC]"), 20};
    };
  // clang-format on

  explicit BaselineAnalyzer(const art::EDAnalyzer::Table<Config>& config);
  void beginRun(art::Run const& run) override;
  void analyze(art::Event const& event) override;
  void endJob() override;
  void FitHistograms();
  void FillSubrunGraphs();
  void WriteThresholds(size_t iCombo);
  void WriteReport();

private:
  // One combination of threshold parameters (one set of outputs each)
  struct ThresholdCombo {
    int offset;
    int offsetPin;
    double nsigma;
  };
  std::string comboTag(size_t iCombo) const;
  std::string comboLabel(size_t iCombo) const;
  static std::string addSuffix(const std::string& fname, const std::string& suffix);

  std::string caloDigiTag_;
  int verbosity_;
  bool writeCSV_;
  std::string CSVfilename_;
  bool writePDF_;
  std::string PDFfilename_;
  std::vector<int> thresholdOffsets_;
  std::vector<int> thresholdOffsetPins_;
  std::vector<double> Nsigmas_;
  std::vector<ThresholdCombo> combos_;   // all combinations, the first one uses the first values
  std::vector<std::string> CSVfilenames_; // one per combination
  double hotStdDev_;
  double coldStdDev_;
  int maxSamples_;
  int minSubrunSamples_;
  int maxPeakToPeak_;
  double minimumStdDev_;
  int minimumThreshold_;

  mu2e::ProditionsHandle<mu2e::CalDAQMap> _calodaqconds_h;

  double xmin;
  double xmax;
  bool titlesSet;

  std::map<int, std::map<int, mu2e::CaloConst::detType>> channelType;
  std::map<int, std::map<int, TH1D*>> h1_baseline_map;

  // Fit results of successfully fitted channels: [board][channel] -> (mean, sigma)
  std::map<int, std::map<int, std::pair<double, double>>> fit_results;
  // Baselines of all channels of the boards with data (defaults for unfitted channels)
  std::map<int, std::map<int, float>> all_baselines;
  // Sample sums per subrun: [subrun][board*100 + channel] -> (sum, count)
  std::map<unsigned, std::map<int, std::pair<double, long>>> subrun_sums;

  TH2D* h2_baselines;
  TH1D* h1_means;
  TH1D* h1_sigmas;
  std::vector<TH1D*> h1_threshold; // one per offset
  TGraphErrors* g_baselines;
  TGraph* g_fitsigmas;
  TGraph* g_stddev;
  TGraph* g_stddev_pin;
  TGraph* g_stddev_lyso;
  TGraph* g_stddev_empty;
  std::vector<TGraph*> g_thresholds; // one per offset
  std::map<int, TGraph*> g_subrun_baselines; // [board*100 + channel], average baseline vs subrun
  TGraph* g_subrun_average;                  // average over all channels vs subrun
  mu2e::THMu2eCaloDisk* h2_disk0;
  mu2e::THMu2eCaloDisk* h2_disk1;
};
} // namespace mu2e

mu2e::BaselineAnalyzer::BaselineAnalyzer(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}, caloDigiTag_(config().caloDigiTag()), verbosity_(config().verbosity()),
    writeCSV_(config().writeCSV()), CSVfilename_(config().CSVfilename()),
    writePDF_(config().writePDF()), PDFfilename_(config().PDFfilename()),
    thresholdOffsets_(config().thresholdOffset()),
    thresholdOffsetPins_(config().thresholdOffsetPin()), Nsigmas_(config().Nsigma()),
    hotStdDev_(config().hotStdDev()),
    coldStdDev_(config().coldStdDev()),
    maxSamples_(config().maxSamples()),
    minSubrunSamples_(config().minSubrunSamples()),
    maxPeakToPeak_(config().maxPeakToPeak()),
    minimumStdDev_(config().minimumStdDev()),
    minimumThreshold_(config().minimumThreshold()) {
  if (thresholdOffsets_.empty() || thresholdOffsetPins_.empty() || Nsigmas_.empty()) {
    throw art::Exception(art::errors::Configuration)
        << "BaselineAnalyzer: thresholdOffset, thresholdOffsetPin and Nsigma must contain at least "
           "one value\n";
  }
  for (int offset : thresholdOffsets_) {
    for (int offsetPin : thresholdOffsetPins_) {
      for (double nsigma : Nsigmas_) {
        combos_.push_back({offset, offsetPin, nsigma});
      }
    }
  }
  art::ServiceHandle<art::TFileService> tfs;

  xmin = 2048 - 150;
  xmax = 2048 + 300;
  h2_baselines = tfs->make<TH2D>("h2_baselines", "Baselines;BoardID*100 + ChannelID;Baseline [ADC]",
                                 CaloConst::_nDIRAC * 100, 0, CaloConst::_nDIRAC * 100, xmax - xmin,
                                 xmin, xmax);
  h1_means = tfs->make<TH1D>("h1_means", "All channels baselines;ADC", int(0.5 * (xmax - xmin)),
                             xmin, xmax);
  h1_sigmas = tfs->make<TH1D>("h1_sigmas", "All channels baseline sigmas;ADC", 100, 0, 10);
  for (size_t i = 0; i < combos_.size(); i++) {
    h1_threshold.push_back(tfs->make<TH1D>(
        ("h1_threshold" + comboTag(i)).c_str(),
        Form("All channel thresholds (%s);ADC", comboLabel(i).c_str()),
        int(0.5 * (xmax - xmin)), xmin, xmax));
    TGraph* g = tfs->makeAndRegister<TGraph>(
        "g_thresholds" + comboTag(i),
        Form("Thresholds (%s);BoardID*100 + ChannelID;ADC", comboLabel(i).c_str()));
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kRed);
    g_thresholds.push_back(g);
  }
  g_baselines = tfs->makeAndRegister<TGraphErrors>(
      "g_baselines", "Baselines;BoardID*100 + ChannelID;Baseline [ADC]");
  g_baselines->SetMarkerStyle(20);
  g_fitsigmas = tfs->makeAndRegister<TGraph>("g_fitsigmas",
                                             "Gaussian fit sigmas;BoardID*100 + ChannelID;ADC");
  g_fitsigmas->SetMarkerStyle(20);
  g_stddev = tfs->makeAndRegister<TGraph>("g_stddev", "Channel StdDev;BoardID*100 + ChannelID;ADC");
  g_stddev->SetMarkerStyle(20);
  g_stddev_pin =
      tfs->makeAndRegister<TGraph>("g_stddev_pin", "[pin-diode];BoardID*100 + ChannelID;ADC");
  g_stddev_pin->SetMarkerStyle(20);
  g_stddev_pin->SetMarkerColor(kRed);
  g_stddev_lyso =
      tfs->makeAndRegister<TGraph>("g_stddev_lyso", "[LYSO];BoardID*100 + ChannelID;ADC");
  g_stddev_lyso->SetMarkerStyle(20);
  g_stddev_lyso->SetMarkerColor(kGreen);
  g_stddev_empty =
      tfs->makeAndRegister<TGraph>("g_stddev_empty", "[empty];BoardID*100 + ChannelID;ADC");
  g_stddev_empty->SetMarkerStyle(20);
  g_stddev_empty->SetMarkerColor(kGray);
  h2_disk0 =
      tfs->makeAndRegister<mu2e::THMu2eCaloDisk>("h2_disk0", "Disk 0", "h2_disk0", "Disk 0", 0);
  h2_disk1 =
      tfs->makeAndRegister<mu2e::THMu2eCaloDisk>("h2_disk1", "Disk 1", "h2_disk1", "Disk 1", 1);
  h2_disk0->SetCombineMode(mu2e::ECombineMode::kAverage);
  h2_disk1->SetCombineMode(mu2e::ECombineMode::kAverage);
  h2_disk0->SetDrawOption("COLZL");
  h2_disk1->SetDrawOption("COLZL");

  auto baseBoardDir = tfs->mkdir("boards");
  for (int boardID = 0; boardID < CaloConst::_nDIRAC; boardID++) {
    TString dirname(Form("Board%03d", boardID));
    auto boardDir = baseBoardDir.mkdir(dirname.Data());
    for (int chanID = 0; chanID < CaloConst::_nChPerDIRAC; chanID++) {
      TString hname = Form("h1_baseline_b%03d_c%02d", boardID, chanID);
      TString htitle = Form("Baseline of Board %03d channel %02d", boardID, chanID);
      h1_baseline_map[boardID][chanID] =
          boardDir.make<TH1D>(hname, htitle, xmax - xmin, xmin, xmax);
    }
  }
  titlesSet = false;
}

// Suffix distinguishing the outputs of the different combinations: only the parameters with
// several values appear (empty if there is only one combination)
std::string mu2e::BaselineAnalyzer::comboTag(size_t iCombo) const {
  const ThresholdCombo& combo = combos_[iCombo];
  std::string tag;
  if (thresholdOffsets_.size() > 1)
    tag += "_offset" + std::to_string(combo.offset);
  if (thresholdOffsetPins_.size() > 1)
    tag += "_offsetPin" + std::to_string(combo.offsetPin);
  if (Nsigmas_.size() > 1)
    tag += Form("_Nsigma%g", combo.nsigma);
  return tag;
}

// Human readable description of a combination, for plot titles
std::string mu2e::BaselineAnalyzer::comboLabel(size_t iCombo) const {
  const ThresholdCombo& combo = combos_[iCombo];
  std::string label = Form("offset %d, pin offset %d", combo.offset, combo.offsetPin);
  if (combo.nsigma != 0)
    label += Form(", %g#sigma", combo.nsigma);
  return label;
}

// Insert suffix before the file extension: "dir/name.ext" -> "dir/name<suffix>.ext"
std::string mu2e::BaselineAnalyzer::addSuffix(const std::string& fname,
                                              const std::string& suffix) {
  size_t slash = fname.find_last_of('/');
  size_t dot = fname.find_last_of('.');
  if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
    return fname + suffix;
  return fname.substr(0, dot) + suffix + fname.substr(dot);
}

void mu2e::BaselineAnalyzer::beginRun(art::Run const& run) {
  if (!CSVfilenames_.empty())
    return; // file names already set by a previous run
  if (PDFfilename_ == "") {
    PDFfilename_ = "baselines_run" + std::to_string(run.id().run()) + ".pdf";
  }
  for (size_t i = 0; i < combos_.size(); i++) {
    std::string defaultName = "thresholds_run" + std::to_string(run.id().run()) + "_offset" +
                              std::to_string(combos_[i].offset) + "_offsetPin" +
                              std::to_string(combos_[i].offsetPin);
    if (combos_[i].nsigma != 0) {
      defaultName += Form("_Nsigma%g", combos_[i].nsigma);
    }
    defaultName += Form("_minRMS%g_minTh%d", minimumStdDev_, minimumThreshold_);
    CSVfilenames_.push_back(CSVfilename_ == "" ? defaultName + ".csv"
                                               : addSuffix(CSVfilename_, comboTag(i)));
  }
}

void mu2e::BaselineAnalyzer::analyze(art::Event const& event) {
  const auto& caloDigis = *event.getValidHandle(consumes<mu2e::CaloDigiCollection>(caloDigiTag_));
  art::ServiceHandle<art::TFileService> tfs;

  mu2e::CalDAQMap const& calodaqconds = _calodaqconds_h.get(event.id());
  if (!titlesSet) { // Need to do this here in order to have calodaqconds
    for (int sipmid = 0; sipmid < CaloConst::_nChannel; sipmid++) {
      mu2e::CaloSiPMId SiPMID_(sipmid);
      if (!SiPMID_.isValid())
        continue;
      mu2e::CaloRawSiPMId rawId = calodaqconds.rawId(SiPMID_);
      int boardID = rawId.dirac();
      int chanID = rawId.ROCchannel();
      channelType[boardID][chanID] = SiPMID_.detType();
      TString typeName = mu2e::CaloConst::detTypeName(channelType[boardID][chanID]);
      TString htitle =
          Form("Baseline of Board %03d channel %02d [%s]", boardID, chanID, typeName.Data());
      h1_baseline_map[boardID][chanID]->SetTitle(htitle);
    }
    titlesSet = true;
  }

  auto& subrun_sum = subrun_sums[event.subRun()];

  // Loop over the calo digis of this event
  for (uint ihit = 0; ihit < caloDigis.size(); ihit++) {
    int SiPMID = caloDigis[ihit].SiPMID();
    std::vector<int> waveform = caloDigis[ihit].waveform();

    mu2e::CaloSiPMId SiPMID_(SiPMID);
    if (!SiPMID_.isValid())
      continue;
    mu2e::CaloRawSiPMId rawId = calodaqconds.rawId(SiPMID_);
    int boardID = rawId.dirac();
    int chanID = rawId.ROCchannel();

    // Discard empty waveforms and cosmics candidates
    if (waveform.empty()) continue;
    if (maxSamples_ < 0){
      auto [wf_min, wf_max] = std::minmax_element(waveform.begin(), waveform.end());
      if (*wf_max - *wf_min >= maxPeakToPeak_) continue;
    }

    // Fill hist
    double sum = 0;
    uint nSamples = (maxSamples_ > 0) ? std::min((int)waveform.size(), maxSamples_) : waveform.size();
    for (uint i = 0; i < nSamples; i++) {
      h2_baselines->Fill(boardID * 100 + chanID, waveform[i]);
      h1_baseline_map[boardID][chanID]->Fill(waveform[i]);
      sum += waveform[i];
    }
    auto& chan_sum = subrun_sum[boardID * 100 + chanID];
    chan_sum.first += sum;
    chan_sum.second += nSamples;
  }
}

void mu2e::BaselineAnalyzer::endJob() {
  // Remove empty hists
  for (int boardID = 0; boardID < CaloConst::_nDIRAC; boardID++) {
    for (int chanID = 0; chanID < CaloConst::_nChPerDIRAC; chanID++) {
      if (h1_baseline_map[boardID][chanID]->GetEntries() == 0) {
        delete h1_baseline_map[boardID][chanID];
        h1_baseline_map[boardID].erase(chanID);
      }
    }
    if (h1_baseline_map[boardID].empty()) {
      h1_baseline_map.erase(boardID);
    }
  }

  FitHistograms();
  FillSubrunGraphs();
  for (size_t i = 0; i < combos_.size(); i++) {
    WriteThresholds(i);
  }
  if (writePDF_) {
    WriteReport();
  }
}

void mu2e::BaselineAnalyzer::FitHistograms() {
  // Perform gaussian fits
  int failed_fits = 0;
  int unprecise_fits = 0;
  std::set<std::pair<int, int>> failed_map;
  std::set<std::pair<int, int>> unprecise_map;
  for (auto board_pair : h1_baseline_map) {
    int board = board_pair.first;

    for (auto channel_pair : board_pair.second) {
      int channel = channel_pair.first;

      h1_baseline_map[board][channel]->Scale(1. / h1_baseline_map[board][channel]->Integral());

      double stddev = h1_baseline_map[board][channel]->GetRMS();
      if (verbosity_ > 0) {
        if (stddev > hotStdDev_) {
          std::cout << "Hot channel: Board " << board << " Channel " << channel << " ("
                    << std::setprecision(4) << std::setw(5) << stddev << ") ["
                    << channelType[board][channel] << "]\n";
        }
        if (stddev < coldStdDev_) {
          std::cout << "Cold channel: Board " << board << " Channel " << channel << " ("
                    << std::setprecision(4) << std::setw(5) << stddev << ") ["
                    << channelType[board][channel] << "]\n";
        }
      }

      switch (channelType[board][channel]) {
      case mu2e::CaloConst::detType::CsI:
        g_stddev->AddPoint(board * 100 + channel, stddev);
        break;
      case mu2e::CaloConst::detType::PINDiode:
        g_stddev_pin->AddPoint(board * 100 + channel, stddev);
        break;
      case mu2e::CaloConst::detType::CAPHRI:
        g_stddev_lyso->AddPoint(board * 100 + channel, stddev);
        break;
      case mu2e::CaloConst::detType::Invalid:
        g_stddev_empty->AddPoint(board * 100 + channel, stddev);
        break;
      default:
        std::cout << "Board " << board << " Channel " << channel << " : unknown type "
                  << channelType[board][channel] << "\n";
      }

      if (verbosity_ > 1)
        std::cout << "Fitting board " << board << " channel " << channel << " ... ";
      double guess_max = h1_baseline_map[board][channel]->GetBinContent(
          h1_baseline_map[board][channel]->GetMaximumBin());
      double guess_mean = h1_baseline_map[board][channel]->GetMean();
      double guess_sigma = h1_baseline_map[board][channel]->GetRMS();

      ((TF1*)gROOT->GetFunction("gaus"))->SetParameters(guess_max, guess_mean, guess_sigma);
      TFitResultPtr fit_result =
          h1_baseline_map[board][channel]->Fit("gaus", "SQW", "", xmin, xmax);
      h1_baseline_map[board][channel]->GetFunction("gaus")->SetNpx(1000);
      if (fit_result < 0) {
        if (verbosity_ > 1)
          std::cout << " <---- BAD FIT!!! Status " << fit_result << "\n";
        failed_fits++;
        failed_map.insert(std::make_pair(board, channel));
        continue;
      }
      double fit_mean = fit_result->Parameter(1);
      double fit_sigma = fit_result->Parameter(2);
      double fit_chi2 = fit_result->Chi2();
      double fit_ndf = fit_result->Ndf();

      if (verbosity_ > 1)
        std::cout << fit_mean << " +- " << fit_sigma << " chi2/ndf=" << fit_chi2 << "/" << fit_ndf;
      if (fit_result > 0) {
        if (verbosity_ > 1)
          std::cout << " <---- UNPRECISE! Status " << fit_result;
        unprecise_fits++;
        unprecise_map.insert(std::make_pair(board, channel));
      }
      if (verbosity_ > 1)
        std::cout << "\n";

      h1_means->Fill(fit_mean);
      h1_sigmas->Fill(fit_sigma);
      g_baselines->SetPoint(g_baselines->GetN(), board * 100 + channel, fit_mean);
      g_baselines->SetPointError(g_baselines->GetN() - 1, 0, fit_sigma);
      g_fitsigmas->SetPoint(g_fitsigmas->GetN(), board * 100 + channel, fit_sigma);
      if (board < 80) {
        h2_disk0->FillRaw(board, channel, stddev);
      } else {
        h2_disk1->FillRaw(board, channel, stddev);
      }

      fit_results[board][channel] = std::make_pair(fit_mean, fit_sigma);
      all_baselines[board][channel] = fit_mean;
    }
  }

  // Check for missing channels -> assign default values
  for (auto board_pair : h1_baseline_map) {
    int b = board_pair.first;
    for (int c = 0; c < 20; c++) {
      if (all_baselines[b].find(c) == all_baselines[b].end()) {
        if (verbosity_ > 0) {
          std::cout << "Setting default values for empty Board " << b << " Channel " << c << "\n";
        }
        all_baselines[b][c] = 2048;
      }
    }
  }

  if (verbosity_ > 0) {
    std::cout << "Failed fits: " << failed_fits << "\n";
    for (auto pair : failed_map) {
      std::cout << "Board " << pair.first << " Channel " << pair.second << "\n";
    }

    std::cout << "Unprecise fits: " << unprecise_fits << "\n";
    for (auto pair : unprecise_map) {
      std::cout << "Board " << pair.first << " Channel " << pair.second << "\n";
    }
  }
}

// Average baseline of each channel per subrun (no fit, statistics are too low), and average of
// the channel averages per subrun
void mu2e::BaselineAnalyzer::FillSubrunGraphs() {
  art::ServiceHandle<art::TFileService> tfs;
  auto subrunDir = tfs->mkdir("subruns");
  g_subrun_average = subrunDir.makeAndRegister<TGraph>(
      "g_subrun_average", "Average baseline of all channels;Subrun;Baseline [ADC]");
  g_subrun_average->SetLineColor(kRed);
  g_subrun_average->SetLineWidth(5);
  g_subrun_average->SetLineStyle(kDashed);
  g_subrun_average->SetMarkerStyle(20);
  g_subrun_average->SetMarkerColor(kRed);

  for (const auto& [subrun, chan_sums] : subrun_sums) {
    double sum_means = 0;
    int nchan = 0;
    for (const auto& [key, chan_sum] : chan_sums) {
      if (chan_sum.second < minSubrunSamples_)
        continue;
      double mean = chan_sum.first / chan_sum.second;
      if (g_subrun_baselines.find(key) == g_subrun_baselines.end()) {
        TGraph* g = subrunDir.makeAndRegister<TGraph>(
            Form("g_subrun_b%03d_c%02d", key / 100, key % 100),
            Form("Baseline of Board %03d channel %02d;Subrun;Baseline [ADC]", key / 100, key % 100));
        g->SetLineColor(kBlack);
        g_subrun_baselines[key] = g;
      }
      g_subrun_baselines[key]->AddPoint(subrun, mean);
      sum_means += mean;
      nchan++;
    }
    if (nchan > 0)
      g_subrun_average->AddPoint(subrun, sum_means / nchan);
  }
}

// Compute the thresholds for one combination of parameters, fill the threshold plots and write
// CSV output
void mu2e::BaselineAnalyzer::WriteThresholds(size_t iCombo) {
  const int offset = combos_[iCombo].offset;
  const int offsetPin = combos_[iCombo].offsetPin;
  const double nsigma = combos_[iCombo].nsigma;
  const std::string& CSVfilename = CSVfilenames_[iCombo];

  bool writeCSV = writeCSV_;
  std::ofstream outputCSV;
  if (writeCSV) {
    outputCSV.open(CSVfilename);
    if (!outputCSV.is_open()) {
      std::cout << "Warning! Can't open file " << CSVfilename << "\n";
      writeCSV = false;
    }
  }

  std::map<int, std::map<int, int>> all_thresholds;
  for (const auto& board_pair : fit_results) {
    int board = board_pair.first;
    for (const auto& channel_pair : board_pair.second) {
      int channel = channel_pair.first;
      double fit_mean = channel_pair.second.first;

      bool pindiode = (channelType[board][channel] == mu2e::CaloConst::detType::PINDiode);
      double rms = std::max(h1_baseline_map[board][channel]->GetRMS(), minimumStdDev_);

      double above = std::max(nsigma * rms + (pindiode ? offsetPin : offset), double(minimumThreshold_));
      int threshold = int(round(fit_mean + above));

      h1_threshold[iCombo]->Fill(threshold);
      g_thresholds[iCombo]->SetPoint(g_thresholds[iCombo]->GetN(), board * 100 + channel,
                                     threshold);
      all_thresholds[board][channel] = threshold;
    }
  }

  // Missing channels -> default thresholds
  for (const auto& board_pair : all_baselines) {
    int b = board_pair.first;
    for (int c = 0; c < 20; c++) {
      if (all_thresholds[b].find(c) == all_thresholds[b].end()) {
        bool pindiode = (channelType[b][c] == mu2e::CaloConst::detType::PINDiode);
        all_thresholds[b][c] = all_baselines[b][c] + std::max(nsigma * minimumStdDev_ + (pindiode ? offsetPin : offset), double(minimumThreshold_));
      }
    }
  }

  if (writeCSV) {
    for (const auto& board_pair : all_baselines) {
      int b = board_pair.first;
      std::stringstream ss;
      ss << b << ",[[";
      for (int c = 0; c < 20; c++) {
        if (c)
          ss << ",";
        ss << all_baselines[b][c];
      }
      ss << "]],[[";
      for (int c = 0; c < 20; c++) {
        if (c)
          ss << ",";
        ss << all_thresholds[b][c];
      }
      ss << "]],\"\",\"\",\"\"";
      std::string csv_line = ss.str();
      outputCSV << csv_line << "\n";
    }
    outputCSV.close();
    std::cout << "wrote file " << CSVfilename << "\n";
  }
}

void mu2e::BaselineAnalyzer::WriteReport() {
  TCanvas* can = new TCanvas("pdf_canvas", "Baseline Report", 2000, 2000);
  can->SetBatch(kTRUE);

  // Silence ROOT Info messages ("Current canvas added to pdf file ..."). art installs its own
  // ROOT error handler which ignores gErrorIgnoreLevel, so filter them out before it sees them.
  static ErrorHandlerFunc_t artErrorHandler = nullptr;
  artErrorHandler = SetErrorHandler([](int level, Bool_t abort, const char* location, const char* msg) {
    if (level >= kWarning && artErrorHandler)
      artErrorHandler(level, abort, location, msg);
  });

  TString pdfname(PDFfilename_.c_str());
  can->SaveAs(pdfname + "[");

  // Summary pages
  can->Clear();
  can->Divide(2, 2);
  can->cd(1);
  h2_baselines->Draw("COLZ");
  gPad->SetLogz();
  can->cd(2);
  h1_means->Draw();
  can->cd(3);
  h1_sigmas->Draw();
  can->cd(4);
  // Only the thresholds of the first combination (first value of each list) are plotted
  h1_threshold[0]->Draw();
  can->SaveAs(pdfname);

  TText* text0 = new TText(0, 50, "DISK 0");
  TText* text1 = new TText(0, 50, "DISK 1");
  TText* textL = new TText(0, -50, "LEFT");
  TText* textR = new TText(0, -50, "RIGHT");
  text0->SetTextAlign(22);
  text1->SetTextAlign(22);
  textL->SetTextAlign(22);
  textR->SetTextAlign(22);

  h2_disk0->GetZaxis()->SetRangeUser(coldStdDev_, hotStdDev_);
  h2_disk1->GetZaxis()->SetRangeUser(coldStdDev_, hotStdDev_);

  gStyle->SetPalette(kViridis);
  can->Clear();
  can->Divide(2, 2);
  can->cd(1);
  h2_disk0->SetCombineMode(mu2e::ECombineMode::kLeft);
  h2_disk0->DrawCopy("colz l");
  text0->Draw("same");
  textL->Draw("same");
  can->cd(2);
  h2_disk0->SetCombineMode(mu2e::ECombineMode::kRight);
  h2_disk0->DrawCopy("colz l");
  text0->Draw("same");
  textR->Draw("same");
  can->cd(3);
  h2_disk1->SetCombineMode(mu2e::ECombineMode::kLeft);
  h2_disk1->DrawCopy("colz l");
  text1->Draw("same");
  textL->Draw("same");
  can->cd(4);
  h2_disk1->SetCombineMode(mu2e::ECombineMode::kRight);
  h2_disk1->DrawCopy("colz l");
  text1->Draw("same");
  textR->Draw("same");
  can->SaveAs(pdfname);

  //Back to kAverage for root file output
  h2_disk0->SetCombineMode(mu2e::ECombineMode::kAverage);
  h2_disk1->SetCombineMode(mu2e::ECombineMode::kAverage);

  // All following pages are 4:3 (only the first two are square)
  can->SetCanvasSize(2000, 1500);

  // Stddev
  can->Clear();
  g_stddev->GetXaxis()->SetLimits(0, CaloConst::_nDIRAC * 100);
  double xmin = g_stddev->GetXaxis()->GetXmin();
  double xmax = g_stddev->GetXaxis()->GetXmax();
  double ymax = g_stddev->GetYaxis()->GetXmax();

  if (g_stddev_pin->GetYaxis()->GetXmax() > ymax) {
    ymax = g_stddev_pin->GetYaxis()->GetXmax();
  }
  if (g_stddev_lyso->GetYaxis()->GetXmax() > ymax) {
    ymax = g_stddev_lyso->GetYaxis()->GetXmax();
  }
  if (hotStdDev_ > ymax) {
    ymax = hotStdDev_;
  }
  g_stddev->GetYaxis()->SetRangeUser(0, ymax + 0.1);

  g_stddev->Draw("AP");
  if (g_stddev_pin->GetN() > 0)
    g_stddev_pin->Draw("P");
  if (g_stddev_lyso->GetN() > 0)
    g_stddev_lyso->Draw("P");
  if (g_stddev_empty->GetN() > 0)
    g_stddev_empty->Draw("P");
  TLegend* leg = new TLegend(0.65, 0.91, 0.9, 0.98);
  leg->AddEntry(g_stddev_pin, g_stddev_pin->GetTitle(), "P");
  leg->AddEntry(g_stddev_lyso, g_stddev_lyso->GetTitle(), "P");
  leg->AddEntry(g_stddev_empty, g_stddev_empty->GetTitle(), "P");
  leg->Draw();

  TLine* l_cold = new TLine(xmin, coldStdDev_, xmax, coldStdDev_);
  TLine* l_hot = new TLine(xmin, hotStdDev_, xmax, hotStdDev_);
  l_cold->SetLineColor(kBlue);
  l_hot->SetLineColor(kRed);
  l_cold->SetLineWidth(2);
  l_hot->SetLineWidth(2);
  l_cold->SetLineStyle(kDashed);
  l_hot->SetLineStyle(kDashed);
  l_cold->Draw("SAME");
  l_hot->Draw("SAME");
  TLatex* text_cold = new TLatex(xmax, coldStdDev_, "#downarrow COLD");
  TLatex* text_hot = new TLatex(xmax, hotStdDev_, "#uparrow HOT");
  text_cold->SetTextSize(0.04);
  text_hot->SetTextSize(0.04);
  text_cold->SetTextAlign(12);
  text_hot->SetTextAlign(12);
  text_cold->SetTextColor(kBlue);
  text_hot->SetTextColor(kRed);
  text_cold->Draw("SAME");
  text_hot->Draw("SAME");
  can->SaveAs(pdfname);
  delete leg;

  // Baselines vs subrun (top) and deviations from the first subrun of each channel (bottom).
  // Only the channels deviating the most are drawn, together with the average of all channels,
  // and listed in the right margin.
  if (g_subrun_average->GetN() > 0) {
    // Copies shifted by the first point of each graph
    auto shifted = [](const TGraph* g) {
      TGraph* gs = (TGraph*)g->Clone();
      gs->SetBit(kCanDelete);
      for (int i = 0; i < gs->GetN(); i++) {
        gs->GetY()[i] -= g->GetY()[0];
      }
      return gs;
    };

    // Largest deviation from the first subrun of each channel: (|deviation|, deviation, key)
    std::vector<std::tuple<double, double, int>> deviations;
    std::map<int, TGraph*> g_shifted;
    for (const auto& [key, g] : g_subrun_baselines) {
      g_shifted[key] = shifted(g);
      double maxdev = 0;
      for (int i = 0; i < g->GetN(); i++) {
        double dev = g_shifted[key]->GetY()[i];
        if (std::abs(dev) > std::abs(maxdev))
          maxdev = dev;
      }
      deviations.emplace_back(std::abs(maxdev), maxdev, key);
    }
    std::sort(deviations.begin(), deviations.end(),
              [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });
    const size_t nDeviating = std::min<size_t>(10, deviations.size());
    const int dev_colors[10] = {kBlue,     kGreen + 2,  kMagenta + 1, kOrange + 7, kCyan + 1,
                                kGray + 2, kYellow + 2, kViolet + 1,  kTeal + 3,   kOrange + 3};
    TGraph* g_avg_shifted = shifted(g_subrun_average);

    const double* sr = g_subrun_average->GetX();
    double sr_lo = sr[0] - 0.5;
    double sr_hi = sr[g_subrun_average->GetN() - 1] + 0.5;

    // Fixed y range, expanded only if a drawn line goes outside
    struct SubrunPad {
      const char* title;
      std::map<int, TGraph*>* graphs;
      TGraph* average;
      double ymin;
      double ymax;
    };
    std::vector<SubrunPad> subrunPads = {
        {"Average baseline per subrun;Subrun;Baseline [ADC]", &g_subrun_baselines,
         g_subrun_average, 2000, 2100},
        {"Baseline deviation from first subrun;Subrun;#Delta baseline [ADC]", &g_shifted,
         g_avg_shifted, -15, 15}};

    can->Clear();
    can->Divide(1, 2);
    for (size_t ipad = 0; ipad < subrunPads.size(); ipad++) {
      const SubrunPad& sp = subrunPads[ipad];
      can->cd(ipad + 1);
      gPad->SetGrid();
      gPad->SetRightMargin(0.25);

      std::vector<TGraph*> drawn = {sp.average};
      for (size_t i = 0; i < nDeviating; i++) {
        drawn.push_back(sp.graphs->at(std::get<2>(deviations[i])));
      }
      double ylo = sp.ymin;
      double yhi = sp.ymax;
      for (const TGraph* g : drawn) {
        for (int i = 0; i < g->GetN(); i++) {
          ylo = std::min(ylo, g->GetY()[i]);
          yhi = std::max(yhi, g->GetY()[i]);
        }
      }
      double margin = 0.05 * (yhi - ylo);
      if (ylo < sp.ymin)
        ylo -= margin;
      if (yhi > sp.ymax)
        yhi += margin;
      gPad->DrawFrame(sr_lo, ylo, sr_hi, yhi, sp.title);

      TLegend* leg_subrun = new TLegend(0.77, 0.1, 0.995, 0.9);
      leg_subrun->SetBit(kCanDelete);
      leg_subrun->SetBorderSize(0);
      leg_subrun->SetFillStyle(0);
      leg_subrun->SetTextSize(0.035);
      leg_subrun->SetMargin(0.2);
      leg_subrun->SetHeader("Max #Delta [ADC]");

      for (size_t i = 0; i < nDeviating; i++) {
        const auto& [absdev, dev, key] = deviations[i];
        TGraph* g = (TGraph*)sp.graphs->at(key)->Clone();
        g->SetBit(kCanDelete);
        g->SetLineColor(dev_colors[i]);
        g->SetMarkerColor(dev_colors[i]);
        g->SetMarkerStyle(20);
        g->SetLineWidth(2);
        g->Draw(g->GetN() > 1 ? "L" : "P");
        leg_subrun->AddEntry(g, Form("B%03d C%02d (%+.1f)", key / 100, key % 100, dev), "L");
      }
      sp.average->Draw(sp.average->GetN() > 1 ? "L" : "P");
      leg_subrun->AddEntry(sp.average, "All channels average", "L");
      leg_subrun->Draw();
    }
    can->SaveAs(pdfname);
  }

  // Baselines map with fitted baselines (red points) and thresholds of the first combination (red lines),
  // one panel per board with data, 8 panels per page
  {
    const int nPanels = 8;
    TLine line;
    line.SetLineColor(kRed);
    line.SetLineWidth(1);
    const double* thr_x = g_thresholds[0]->GetX();
    const double* thr_y = g_thresholds[0]->GetY();
    int panel = 0;
    for (const auto& board_pair : h1_baseline_map) {
      int board = board_pair.first;
      if (panel % nPanels == 0) {
        can->Clear();
        can->Divide(4, 2);
      }
      can->cd(panel % nPanels + 1);
      gPad->SetLogz();
      TH1* h = h2_baselines->DrawCopy("COL");
      h->SetTitle(Form("Board %d", board));
      h->GetXaxis()->SetRangeUser(board * 100, board * 100 + 20);
      h->GetYaxis()->SetRangeUser(1950, 2200);
      TGraphErrors* g_baselines_draw = (TGraphErrors*)g_baselines->DrawClone("P X");
      g_baselines_draw->SetMarkerColor(kRed);
      g_baselines_draw->SetMarkerSize(0.8);
      for (int i = 0; i < g_baselines_draw->GetN(); i++) {
        g_baselines_draw->GetX()[i] += 0.5; // center of the channel bin
      }
      for (int i = 0; i < g_thresholds[0]->GetN(); i++) {
        if (thr_x[i] >= board * 100 && thr_x[i] < board * 100 + 20) {
          line.DrawLine(thr_x[i], thr_y[i], thr_x[i] + 1, thr_y[i]); // one bin wide
        }
      }
      // Legend in the empty band at the top of the plot (deleted with the pad on Clear)
      TLegend* leg_map = new TLegend(0.12, 0.80, 0.88, 0.89);
      leg_map->SetBit(kCanDelete);
      leg_map->SetBorderSize(0);
      leg_map->SetFillStyle(0);
      leg_map->SetTextSize(0.035);
      leg_map->AddEntry(g_baselines_draw, "Fitted baseline", "P");
      leg_map->AddEntry(&line, Form("Threshold (%s)", comboLabel(0).c_str()), "L");
      leg_map->Draw();
      panel++;
      if (panel % nPanels == 0) {
        can->SaveAs(pdfname);
      }
    }
    if (panel % nPanels != 0) {
      can->SaveAs(pdfname);
    }
  }

  // Hottest channels (largest StdDev) overlaid
  const size_t nHottest = 10;
  std::vector<std::tuple<double, int, int>> channel_stddevs; // (stddev, board, channel)
  for (const auto& board_pair : h1_baseline_map) {
    for (const auto& channel_pair : board_pair.second) {
      channel_stddevs.emplace_back(channel_pair.second->GetRMS(), board_pair.first,
                                   channel_pair.first);
    }
  }
  std::sort(channel_stddevs.begin(), channel_stddevs.end(),
            [](const auto& a, const auto& b) { return std::get<0>(a) > std::get<0>(b); });
  if (channel_stddevs.size() > nHottest) {
    channel_stddevs.resize(nHottest);
  }
  if (!channel_stddevs.empty()) {
    const int colors[nHottest] = {kRed,     kBlue,      kGreen + 2, kMagenta, kCyan + 1,
                                  kOrange + 7, kViolet + 1, kAzure + 7, kPink + 9, kGray + 2};
    can->Clear();
    can->SetLogy();
    TLegend* leg_hot = new TLegend(0.65, 0.55, 0.9, 0.9);
    leg_hot->SetHeader(Form("%zu hottest channels", channel_stddevs.size()));
    double hot_max = 0;
    for (const auto& [stddev, board, channel] : channel_stddevs) {
      hot_max = std::max(hot_max, h1_baseline_map[board][channel]->GetMaximum());
    }
    for (size_t i = 0; i < channel_stddevs.size(); i++) {
      const auto& [stddev, board, channel] = channel_stddevs[i];
      // Draw copies so the per-board pages below keep the default style
      TH1* h = h1_baseline_map[board][channel]->DrawCopy(i == 0 ? "HIST" : "HIST SAME");
      h->SetTitle("Hottest channels;ADC");
      h->GetXaxis()->SetRangeUser(1950, 2150);
      h->SetMaximum(2 * hot_max);
      h->SetMinimum(1e-6);
      h->SetLineColor(colors[i]);
      h->SetLineWidth(2);
      leg_hot->AddEntry(h, Form("B%03d C%02d (%.2f) [%s]", board, channel, stddev,
                                mu2e::CaloConst::detTypeName(channelType[board][channel]).c_str()),
                        "L");
    }
    leg_hot->Draw();
    can->SaveAs(pdfname);
    delete leg_hot;
    can->SetLogy(0);
  }

  // Per-board channel histograms, grouped into pages of 20 plots each
  for (const auto& board_pair : h1_baseline_map) {
    int board = board_pair.first;
    can->Clear();
    can->Divide(5, 4);
    for (int channel = 0; channel < 20; channel++) {
      can->cd(channel + 1);
      if (h1_baseline_map[board].find(channel) != h1_baseline_map[board].end()) {
        h1_baseline_map[board][channel]->GetXaxis()->SetRangeUser(1950, 2150);
        h1_baseline_map[board][channel]->Draw("HIST");
        gPad->SetLogy();
      }
    }
    can->SaveAs(pdfname);
  }

  can->SaveAs(pdfname + "]");
  delete can;

  SetErrorHandler(artErrorHandler);
}

DEFINE_ART_MODULE(mu2e::BaselineAnalyzer)
