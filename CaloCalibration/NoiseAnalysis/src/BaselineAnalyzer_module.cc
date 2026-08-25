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
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TEllipse.h"
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
      fhicl::Atom<bool> writeTXT {fhicl::Name("writeTXT" ) , fhicl::Comment("Write per-board text files with thresholds"), false};
      fhicl::Atom<std::string> TXTfoldername {fhicl::Name("TXTfoldername" ) , fhicl::Comment("Folder to write thresholds into"), ""};
      fhicl::Atom<bool> writeCSV {fhicl::Name("writeCSV" ) , fhicl::Comment("Write CSV file with thresholds"), false};
      fhicl::Atom<std::string> CSVfilename {fhicl::Name("CSVfilename" ) , fhicl::Comment("CSV file to write thresholds into"), ""};
      fhicl::Atom<bool> writePDF {fhicl::Name("writePDF" ) , fhicl::Comment("Write PDF report"), false};
      fhicl::Atom<std::string> PDFfilename {fhicl::Name("PDFfilename" ) , fhicl::Comment("PDF report file name"), ""};
      fhicl::Atom<int> thresholdOffset {fhicl::Name("thresholdOffset" ) , fhicl::Comment("Offset with respect to gaussian mean"), 100};
      fhicl::Atom<int> thresholdOffsetPin {fhicl::Name("thresholdOffsetPin" ) , fhicl::Comment("Offset with respect to gaussian mean (pin diodes)"), 100};
      fhicl::Atom<int> hotStdDev {fhicl::Name("hotStdDev" ) , fhicl::Comment("StdDev limit for hot channels"), 6};
      fhicl::Atom<int> coldStdDev {fhicl::Name("coldStdDev" ) , fhicl::Comment("StdDev limit for cold channels"), 3};
    };
  // clang-format on

  explicit BaselineAnalyzer(const art::EDAnalyzer::Table<Config>& config);
  void beginRun(art::Run const& run) override;
  void analyze(art::Event const& event) override;
  void endJob() override;
  void FitHistograms();
  void WriteReport();

private:
  std::string caloDigiTag_;
  int verbosity_;
  bool writeTXT_;
  std::string TXTfoldername_;
  bool writeCSV_;
  std::string CSVfilename_;
  bool writePDF_;
  std::string PDFfilename_;
  int thresholdOffset_;
  int thresholdOffsetPin_;
  double hotStdDev_;
  double coldStdDev_;

  mu2e::ProditionsHandle<mu2e::CalDAQMap> _calodaqconds_h;

  double xmin;
  double xmax;
  bool titlesSet;

  std::map<int, std::map<int, mu2e::CaloConst::detType>> channelType;
  std::map<int, std::map<int, TH1D*>> h1_baseline_map;

  TH2D* h2_baselines;
  TH1D* h1_means;
  TH1D* h1_sigmas;
  TH1D* h1_threshold;
  TGraphErrors* g_baselines;
  TGraph* g_fitsigmas;
  TGraph* g_stddev;
  TGraph* g_stddev_pin;
  TGraph* g_stddev_lyso;
  TGraph* g_stddev_empty;
  TGraph* g_thresholds;
  mu2e::THMu2eCaloDisk* h2_disk0;
  mu2e::THMu2eCaloDisk* h2_disk1;
};
} // namespace mu2e

mu2e::BaselineAnalyzer::BaselineAnalyzer(const art::EDAnalyzer::Table<Config>& config) :
    art::EDAnalyzer{config}, caloDigiTag_(config().caloDigiTag()), verbosity_(config().verbosity()),
    writeTXT_(config().writeTXT()), TXTfoldername_(config().TXTfoldername()),
    writeCSV_(config().writeCSV()), CSVfilename_(config().CSVfilename()),
    writePDF_(config().writePDF()), PDFfilename_(config().PDFfilename()),
    thresholdOffset_(config().thresholdOffset()),
    thresholdOffsetPin_(config().thresholdOffsetPin()), hotStdDev_(config().hotStdDev()),
    coldStdDev_(config().coldStdDev()) {
  art::ServiceHandle<art::TFileService> tfs;

  xmin = 2048 - 150;
  xmax = 2048 + 300;
  h2_baselines = tfs->make<TH2D>("h2_baselines", "Baselines;BoardID*100 + ChannelID;Baseline [ADC]",
                                 CaloConst::_nDIRAC*100, 0, CaloConst::_nDIRAC*100, xmax - xmin, xmin, xmax);
  h1_means = tfs->make<TH1D>("h1_means", "All channels baselines;ADC", int(0.5 * (xmax - xmin)),
                             xmin, xmax);
  h1_sigmas = tfs->make<TH1D>("h1_sigmas", "All channels baseline sigmas;ADC", 100, 0, 10);
  h1_threshold = tfs->make<TH1D>("h1_threshold", "All channel thresholds;ADC",
                                 int(0.5 * (xmax - xmin)), xmin, xmax);
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
  g_thresholds =
      tfs->makeAndRegister<TGraph>("g_thresholds", "Thresholds;BoardID*100 + ChannelID;ADC");
  g_thresholds->SetMarkerStyle(20);
  g_thresholds->SetMarkerColor(kRed);
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

void mu2e::BaselineAnalyzer::beginRun(art::Run const& run) {
  if (CSVfilename_ == "") {
    CSVfilename_ = std::string("thresholds_run") + run.id().run() + "_offset" + thresholdOffset_ +
                   "_offsetPin" + thresholdOffsetPin_ + ".csv";
  }
  if (PDFfilename_ == "") {
    PDFfilename_ = std::string("thresholds_run") + run.id().run() + "_offset" + thresholdOffset_ +
                   "_offsetPin" + thresholdOffsetPin_ + ".pdf";
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

    // Fill hist
    for (auto sample : waveform) {
      h2_baselines->Fill(boardID * 100 + chanID, sample);
      h1_baseline_map[boardID][chanID]->Fill(sample);
    }
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
  std::ofstream outputCSV;
  if (writeCSV_) {
    outputCSV.open(CSVfilename_);
    if (!outputCSV.is_open()) {
      std::cout << "Warning! Can't open file " << CSVfilename_ << "\n";
      writeCSV_ = false;
    }
  }
  std::set<int> all_boards;
  std::map<int, std::map<int, float>> all_baselines;
  std::map<int, std::map<int, int>> all_thresholds;
  for (auto board_pair : h1_baseline_map) {
    int board = board_pair.first;
    all_boards.insert(board);
    std::ofstream outputBaselineFile;
    TString fname = Form("%s/dirac%03d.baseline", TXTfoldername_.c_str(), board);

    if (writeTXT_) {
      outputBaselineFile.open(fname);
      if (!outputBaselineFile.is_open()) {
        std::cout << "Warning! Can't open file " << fname << "\n";
        writeTXT_ = false;
      }
    }

    for (auto channel_pair : board_pair.second) {
      int channel = channel_pair.first;

      bool pindiode = (channelType[board][channel] == mu2e::CaloConst::detType::PINDiode);

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
      int threshold = int(round(fit_mean + thresholdOffset_));

      // Pin diodes
      if (pindiode) {
        threshold = int(round(fit_mean + thresholdOffsetPin_));
      }

      if (verbosity_ > 0)
        std::cout << fit_mean << " +- " << fit_sigma << " chi2/ndf=" << fit_chi2 << "/" << fit_ndf;
      if (fit_result > 0) {
        if (verbosity_ > 0)
          std::cout << " <---- UNPRECISE! Status " << fit_result;
        unprecise_fits++;
        unprecise_map.insert(std::make_pair(board, channel));
      }
      if (verbosity_ > 0)
        std::cout << "\n";

      h1_means->Fill(fit_mean);
      h1_sigmas->Fill(fit_sigma);
      h1_threshold->Fill(threshold);
      g_baselines->SetPoint(g_baselines->GetN(), board * 100 + channel, fit_mean);
      g_baselines->SetPointError(g_baselines->GetN() - 1, 0, fit_sigma);
      g_fitsigmas->SetPoint(g_fitsigmas->GetN(), board * 100 + channel, fit_sigma);
      g_thresholds->SetPoint(g_thresholds->GetN(), board * 100 + channel, threshold);
      if (board < 80) {
        h2_disk0->FillRaw(board, channel, stddev);
      } else {
        h2_disk1->FillRaw(board, channel, stddev);
      }

      all_baselines[board][channel] = fit_mean;
      all_thresholds[board][channel] = threshold;

      TString lineout = Form("%d\t%.2f\t%.2f\t%d", channel, fit_mean, fit_sigma, threshold);
      if (writeTXT_) {
        outputBaselineFile << lineout << "\n";
      }
    }
    if (writeTXT_) {
      outputBaselineFile.close();
    }
  }

  // Check for missing channels -> assign default values
  for (auto b : all_boards) {
    for (int c = 0; c < 20; c++) {
      if (all_baselines[b].find(c) == all_baselines[b].end()) {
        if (verbosity_ > 0) {
          std::cout << "Setting default values for empty Board " << b << " Channel " << c << "\n";
        }
        all_baselines[b][c] = 2048;
        all_thresholds[b][c] = all_baselines[b][c] + thresholdOffset_;
      }
    }
  }

  if (writeCSV_) {
    for (auto b : all_boards) {
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
    std::cout << "wrote file " << CSVfilename_ << "\n";
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

void mu2e::BaselineAnalyzer::WriteReport() {
  TCanvas* can = new TCanvas("pdf_canvas", "Baseline Report", 2000, 1400);
  can->SetBatch(kTRUE);

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
  h1_threshold->Draw();
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
  h2_disk0->Draw("colz l");
  text0->Draw("same");
  textL->Draw("same");
  can->cd(2);
  h2_disk0->SetCombineMode(mu2e::ECombineMode::kRight);
  h2_disk0->Draw("colz l");
  text0->Draw("same");
  textR->Draw("same");
  can->cd(3);
  h2_disk1->SetCombineMode(mu2e::ECombineMode::kLeft);
  h2_disk1->Draw("colz l");
  text1->Draw("same");
  textL->Draw("same");
  can->cd(4);
  h2_disk1->SetCombineMode(mu2e::ECombineMode::kRight);
  h2_disk1->Draw("colz l");
  text1->Draw("same");
  textR->Draw("same");
  can->SaveAs(pdfname);

  // Stddev
  can->Clear();
  g_stddev->GetXaxis()->SetLimits(0, CaloConst::_nDIRAC*100);
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
}

DEFINE_ART_MODULE(mu2e::BaselineAnalyzer)
