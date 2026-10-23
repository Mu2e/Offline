// ROOT styling for DQM plots, shared by the online and offline DQM modules.
// Adapted from mu2e.mplstyle
// Samuel Grant 2025

#ifndef DQMHelpers_inc_DQMStyle_hh
#define DQMHelpers_inc_DQMStyle_hh

#include "TColor.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"

#include <string>

namespace mu2e {

// Styling only: nothing here books, fills or interprets a histogram.
class DQMStyle {
 public:
  // Installs the DQM style process-wide; call once, from whoever owns the drawing.
  static void SetStyle() {
    TStyle* style = new TStyle("Mu2eDQM", "Mu2e DQM styling");

    // White, borderless canvas and pads
    style->SetCanvasColor(kWhite);
    style->SetCanvasBorderMode(0);
    style->SetCanvasDefH(600);
    style->SetCanvasDefW(800);
    style->SetPadColor(kWhite);
    style->SetPadBorderMode(0);

    // Margins
    style->SetPadTopMargin(0.08);
    style->SetPadBottomMargin(0.14);
    style->SetPadLeftMargin(0.14);
    style->SetPadRightMargin(0.06);

    // Helvetica (font 42) for axes, labels, stats
    int font = 42;
    style->SetTextFont(font);
    style->SetLabelFont(font, "xyz");
    style->SetTitleFont(font, "xyz");
    style->SetStatFont(font);

    // Font sizes
    style->SetTextSize(0.045);
    style->SetTitleSize(0.045, "xyz");
    style->SetLabelSize(0.045, "xyz");
    style->SetStatFontSize(0.040);
    style->SetLegendTextSize(0.040);
    style->SetLegendBorderSize(0);
    style->SetLegendFillColor(kWhite);

    // Axis title offsets
    style->SetTitleOffset(1.15, "x");
    style->SetTitleOffset(1.50, "y");
    style->SetTitleOffset(1.15, "z");

    // Title box — borderless, transparent, centred, bold
    style->SetOptTitle(1);
    style->SetTitleBorderSize(0);
    style->SetTitleFillColor(0);
    style->SetTitleStyle(0);
    style->SetTitleFont(62, "");  // Helvetica bold
    style->SetTitleFontSize(0.045);
    style->SetTitleAlign(23);
    style->SetTitleX(0.5);

    // Stat box — entries, mean, RMS, underflow, overflow; frameless
    // Numeric format ksiourmen: 000111110 = 111110
    style->SetOptStat(111110);
    style->SetStatBorderSize(0);
    style->SetStatStyle(0);
    style->SetStatX(0.92);
    style->SetStatY(0.90);
    style->SetStatW(0.22);
    style->SetStatH(0.20);

    // Fills
    style->SetFillColor(kWhite);
    style->SetFillStyle(1001);

    // Default histogram colours (used by ForceStyle)
    style->SetHistLineColor(kBlack);
    style->SetHistFillColor(kGray);
    style->SetHistFillStyle(1001);

    // axes.linewidth: 2
    style->SetFrameLineWidth(2);

    // Ticks on all sides, inward
    style->SetPadTickX(1);
    style->SetPadTickY(1);
    style->SetTickLength(0.03, "x");
    style->SetTickLength(0.03, "y");
    style->SetTickLength(0.03, "z");

    // Minor ticks visible (510 = 10 primary, 5 minor)
    style->SetNdivisions(510, "x");
    style->SetNdivisions(510, "y");
    style->SetNdivisions(510, "z");

    // Line widths
    style->SetLineWidth(2);
    style->SetHistLineWidth(2);

    // Grid
    style->SetGridColor(kGray + 1);
    style->SetGridStyle(3);
    style->SetGridWidth(1);

    // 2D palette
    style->SetPalette(kInvertedDarkBodyRadiator);
    style->SetNumberContours(256);

    gROOT->SetStyle("Mu2eDQM");
    gROOT->ForceStyle();
  }

  static void FormatHist(TH1* hist, const std::string& colour = "black") {
    if (!hist) return;

    hist->SetLineWidth(2);
    hist->SetFillStyle(1001);
    hist->SetLineStyle(1);
    hist->GetYaxis()->SetTitleOffset(1.60);
    hist->GetXaxis()->SetTitleOffset(1.15);

    if (colour == "black") {  // default: black/grey
      hist->SetLineColor(kBlack);
      hist->SetFillColor(kGray);
    } else if (colour == "blue") {
      hist->SetLineColor(kAzure + 2);
      hist->SetFillColor(kAzure - 9);
    } else if (colour == "green") {
      hist->SetLineColor(kGreen + 2);
      hist->SetFillColor(kGreen - 9);
    } else if (colour == "red") {
      hist->SetLineColor(kRed + 2);
      hist->SetFillColor(kRed - 9);
    }
  }

  static void FormatHist2D(TH2* hist) {
    if (!hist) return;

    hist->GetXaxis()->SetTitleOffset(1.15);
    hist->GetYaxis()->SetTitleOffset(1.50);
    hist->GetZaxis()->SetTitleOffset(1.15);
    hist->SetStats(0);
  }

  static void FormatGraph(TGraph* graph, const std::string& colour = "black") {
    if (!graph) return;

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.8);
    graph->SetLineStyle(1);
    graph->GetYaxis()->SetTitleOffset(1.60);
    graph->GetXaxis()->SetTitleOffset(1.15);

    if (colour == "black") {  // default: black/grey
      graph->SetLineColor(kBlack);
      graph->SetFillColor(kGray);
    } else if (colour == "blue") {
      graph->SetLineColor(kAzure + 2);
      graph->SetFillColor(kAzure - 9);
    } else if (colour == "green") {
      graph->SetLineColor(kGreen + 2);
      graph->SetFillColor(kGreen - 9);
    } else if (colour == "red") {
      graph->SetLineColor(kRed + 2);
      graph->SetFillColor(kRed - 9);
    }
  }
};

}  // namespace mu2e

#endif /* DQMHelpers_inc_DQMStyle_hh */
