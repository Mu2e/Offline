// Display functions for the AgnosticHelixFinder diagnostics
// Original authors: Michael MacKenzie, Matthew Stortini

#include "Offline/CalPatRec/inc/AgnosticHelixFinderDiag.hh"

using namespace mu2e;

//-----------------------------------------------------------------------------
// Standard color choices by pdg/momentum
int AgnosticHelixFinderDiag::MCColor(int pdg, double pmc) {
  pdg = std::abs(pdg);
  if(pdg == 0) return kBlack;
  if(pdg == 2212) return kOrange;
  if(pmc > 0.) {
    if(pmc < 20.) return kAzure-9;
    if(pmc < 50.) return kCyan;
    if(pmc < 70.) return kMagenta - 2;
  }
  return kBlue;
}

//-----------------------------------------------------------------------------
// Plot label information for a display pad
void AgnosticHelixFinderDiag::plotLabel(float left_margin, float top_margin,
                                        float right_margin , std::string title) {
  TLatex* label = new TLatex();
  label->SetNDC();
  label->SetTextFont(43);
  label->SetTextSize(24);
  label->SetTextAlign(11);
  label->DrawLatex(left_margin, 1.01 - top_margin,
                   Form("Event: %i:%i:%u",
                        _data->event->run(),
                        _data->event->subRun(),
                        _data->event->event()));
  if(title != "") {
    label->SetTextAlign(31);
    label->DrawLatex(1. - right_margin, 1.01 - top_margin, title.c_str());
  }
  _drawn_objects.push_back(label);
}

//-----------------------------------------------------------------------------
// Initialize data/objects for the display
void AgnosticHelixFinderDiag::initDisplay() {

  // Need the TApplication for the display, where TRint allows command-line usage
  int    tmp_argc(0);
  char** tmp_argv(0);
  if(!gApplication && !_application) {
    _application = new TRint("AgnosticHelixFinderDiag_tool", &tmp_argc, tmp_argv);
    gApplication = _application;
  } else if(gApplication && !_application) _application = gApplication; // use the existing one

  // Ensure the display canvases are available
  if(_display3D && !_c_3d) {
    _c_3d = new TCanvas("c_3d", "3D canvas"  , 900, 900);
  }
  if(!_c_tot) {
    _c_tot = new TCanvas("c_tot", "Total canvas"  , 1260, 900);
  }
  if(!_c_hlx) {
    _c_hlx = new TPad("c_hlx", "Helix canvas"  , 0., 0.3, 0.5, 1.);
    _c_hlx->SetLeftMargin(0.12);
    _c_hlx->SetBottomMargin(0.10);
    _c_hlx->SetTicks(1);
    _c_tot->cd();
    _c_hlx->Draw();
  }
  if(!_c_seg) {
    _c_seg = new TPad("c_seg", "Segment canvas"  , 0., 0.0, 1.0, 0.3);
    _c_seg->SetLeftMargin(0.06);
    _c_seg->SetBottomMargin(0.22);
    _c_seg->SetRightMargin(0.01);
    _c_seg->SetTopMargin(0.08);
    _c_seg->SetTicks(1);
    _c_tot->cd();
    _c_seg->Draw();
  }
  if(!_c_trp) {
    _c_trp = new TPad("c_trp", "Triplet canvas"  , 0.5, 0.3, 1.0, 1.0);
    _c_trp->SetLeftMargin(0.12);
    _c_trp->SetBottomMargin(0.12);
    _c_trp->SetTicks(1);
    _c_tot->cd();
    _c_trp->Draw();
  }

  // clean up the last event
  for(auto o : _drawn_objects) delete o;
  _drawn_objects.clear();
}

//-----------------------------------------------------------------------------
// Draw an x-y circle
void AgnosticHelixFinderDiag::plotCircle(float xc, float yc, float r,
                                         int color, int marker_style,
                                         std::string title,
                                         bool draw_center) {
  // set up graph point for the center point
  if(draw_center) {
    float centerX[1] = {xc};
    float centerY[1] = {yc};
    TGraph* center = new TGraph(1, centerX, centerY);
    center->SetMarkerStyle(marker_style);
    center->SetMarkerColor(color);
    center->Draw("P");
    center->SetName((title + "_center").c_str());
    center->SetTitle((title + " center").c_str());
    _drawn_objects.push_back(center);
  }

  // draw the circle
  TEllipse* circle = new TEllipse(xc, yc, r, r);
  circle->SetLineColor(color);
  circle->SetFillStyle(0);
  circle->Draw("same");
  _drawn_objects.push_back(circle);

}

//-----------------------------------------------------------------------------
// Plot a straw hit
void AgnosticHelixFinderDiag::plotHitXY(const ComboHit& hit, int color) {

  // Get the hit position
  const auto pos = hit.pos();

  // error bar along the wire
  const float wireDirX = hit.uDir().x();
  const float wireDirY = hit.uDir().y();
  const float wireDirErr = hit.wireRes();
  const float xWire1 = pos.x() - wireDirX * wireDirErr;
  const float xWire2 = pos.x() + wireDirX * wireDirErr;
  const float yWire1 = pos.y() - wireDirY * wireDirErr;
  const float yWire2 = pos.y() + wireDirY * wireDirErr;

  // error bar transverse to the wire
  const float perpDirX = wireDirY;
  const float perpDirY = -wireDirX;
  const float perpDirErr = hit.transRes();
  const float xPerp1 = pos.x() - perpDirX * perpDirErr;
  const float xPerp2 = pos.x() + perpDirX * perpDirErr;
  const float yPerp1 = pos.y() - perpDirY * perpDirErr;
  const float yPerp2 = pos.y() + perpDirY * perpDirErr;

  // initialize lines
  auto w_line = new TLine(xWire1, yWire1, xWire2, yWire2);
  w_line->SetLineColor(color);
  w_line->Draw("same");
  auto t_line = new TLine(xPerp1, yPerp1, xPerp2, yPerp2);
  t_line->SetLineColor(color);
  t_line->Draw("same");
  _drawn_objects.push_back(w_line);
  _drawn_objects.push_back(t_line);
}

//-----------------------------------------------------------------------------
// Plot a straw hit associated with a triplet
void AgnosticHelixFinderDiag::plotHitXY(const tripletPoint& point, int color) {
  if(!_data || !_data->chColl) return;
  if(point.hitIndice >= 0) {
    if(point.hitIndice >= int(_data->chColl->size())) return;
    plotHitXY(_data->chColl->at(point.hitIndice), color); // normal straw hit
    return;
  }

  // draw the fake hit as a filled circle
  const auto pos = point.pos;
  TEllipse* circle = new TEllipse(pos.x(), pos.y(), 10., 10.);
  circle->SetLineColor(color);
  circle->SetFillColor(color);
  circle->SetFillStyle(3001);
  circle->Draw("same");
  _drawn_objects.push_back(circle);
}

//-----------------------------------------------------------------------------
// Plot a hit in phi-z
void AgnosticHelixFinderDiag::plotHitPhiZ(double z, double phi, int color, int iline, int marker) {

  // initialize the point
  TGraph* g = new TGraph(1, &z, &phi);
  constexpr int styles[] = {24, 25, 26, 27, 28, 30};
  const static int nstyles = sizeof(styles) / sizeof(*styles);
  if(marker < 0) {
    const int istyle = std::max(iline, 0) % nstyles;
    marker = styles[istyle];
  }
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1);
  g->Draw("PE");
  if(_debugLevel > 1) printf("  %s: Adding point z = %6.1f, phi = %5.2f, line = %i\n", __func__, z, phi, iline);

  _drawn_objects.push_back(g);
}

//-----------------------------------------------------------------------------
// Plot a combo hit in phi-z, for the assumed helix
void AgnosticHelixFinderDiag::plotHitPhiZ(int index, int color, int iline, int marker) {
  if(!_data || !_data->tcHits || !_data->chColl) return;

  // Get the hit
  const auto& hit = _data->tcHits->at(index);

  // Get the hit info
  const float phi = hit.helixPhi; // + 2 * M_PI * line.helixPhiCorrections[i];
  const float phiError = std::sqrt(hit.helixPhiError2);
  const int hit_index = hit.hitIndice;
  float z(0.), zError(0.);

  if(hit_index >= 0) {
    const auto& combo_hit = _data->chColl->at(hit_index);
    z = combo_hit.pos().z();
  } else if(hit_index == HitType::CALOCLUSTER) {
    z = _data->caloPos.z();
  } else if(hit_index == HitType::STOPPINGTARGET) {
    z = _data->caloPos.z();
  }

  // initialize the point
  TGraph* g = new TGraphErrors(1, &z, &phi, &zError, &phiError);
  constexpr int styles[] = {24, 25, 26, 27, 28, 30};
  const static int nstyles = sizeof(styles) / sizeof(*styles);
  if(marker < 0) {
    const int istyle = std::max(iline, 0) % nstyles;
    marker = styles[istyle];
  }
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1);
  g->Draw("PE");
  if(_debugLevel > 1) printf("  %s: Adding point z = %6.1f, phi = %5.2f, line = %i\n", __func__, z, phi, iline);

  _drawn_objects.push_back(g);
}

//-----------------------------------------------------------------------------
// Plot a fit line in phi-z
void AgnosticHelixFinderDiag::plotLinePhiZ(const lineInfo& line, bool resolve, int iline) {

  constexpr int colors[] = {kBlue, kGreen, kRed, kOrange, kPink, kYellow-2};
  const static int ncolors = sizeof(colors) / sizeof(*colors);
  const int color = colors[iline % ncolors];

  // Plot the points
  for(auto index : line.tcHitsIndices) plotHitPhiZ(index, color, iline);

  // Draw the fit line
  const double zmin  = line.zMin;
  const double zmax  = line.zMax;
  // clone the fitter to avoid changing its state
  ::LsqSums2 fitter = line.fitter;
  const double slope = fitter.dydx();
  const int   traj   = (slope < 0.) ? 0 : 1;
  const double phi0  = fitter.y0();;
  const double phi1  = slope * zmin + phi0;
  const double phi2  = slope * zmax + phi0;

  // Add line segments that mod by pi
  double z_start(zmin), phi_start(std::fmod(phi1, 2.*M_PI));
  if(phi_start < 0.) phi_start += 2.*M_PI;
  int max_loops(100), loops(0); // protect against infinite loops
  while(loops < max_loops && z_start < zmax) {
    double z_end = (slope == 0.) ? zmax : std::min(zmax, z_start + (traj*2.*M_PI - phi_start)/slope);
    double phi_end = phi_start + (z_end - z_start)*slope;
    if(resolve) { // use the unchanged values if resolving loop numbers
      z_start = zmin;
      phi_start = phi1;
      z_end = zmax;
      phi_end = phi2;
    }
    TLine* tline = new TLine(z_start, phi_start, z_end, phi_end);
    tline->SetLineColor(color);
    tline->SetLineWidth(1);
    tline->Draw("same");
    _drawn_objects.push_back(tline);

    z_start = z_end;
    phi_start = (1-traj)*2.*M_PI;
    ++loops;
  }
}

//-----------------------------------------------------------------------------
// Plot a helix in x-y
void AgnosticHelixFinderDiag::plotHelixXY(const HelixSeed& hseed, const int index) {
  if(hseed.hits().empty()) return;

  const float r  = hseed.helix().radius(); // reco info
  const float xC = hseed.helix().center().x();
  const float yC = hseed.helix().center().y();
  if(_debugLevel > 2) printf("    Reco circle %i: x = %6.1f, y = %6.1f, r = %6.1f\n", index, xC, yC, r);
  std::string title = (index >= 0) ? Form("Reco_%i", index) : "Reco";
  plotCircle(xC, yC, r, kRed, 2, title); // Reco helix

  // Helix hits
  for(auto& hit : hseed.hits()) plotHitXY(hit, kGreen);

  // check for MC info
  const auto sim = helixSimMatch(hseed);
  const int sim_id = (sim) ? sim->id().asInt() : -1;
  if(sim && _simInfo.find(sim_id) != _simInfo.end()) {
    float xC_mc, yC_mc, r_mc;
    MCCircle(sim_id, xC_mc, yC_mc, r_mc);
    if(_debugLevel > 2) printf("    MC circle   %i: x = %6.1f, y = %6.1f, r = %6.1f\n", index, xC_mc, yC_mc, r_mc);
    title = (index >= 0) ? Form("MC_%i", index) : "MC";
    const int pdg = std::abs(sim->pdgId());
    if(pdg != 2212 || _showProtons) {
      const int mc_color = MCColor(pdg);
      plotCircle(xC_mc, yC_mc, r_mc, mc_color, 2, title); // MC helix estimate
    }
  }
}

//-----------------------------------------------------------------------------
// Plot a helix in phi-z
void AgnosticHelixFinderDiag::plotHelixPhiZ(const HelixSeed& hseed, const int index) {
  if(hseed.hits().empty()) return;

  // Helix hits
  const int iline = std::max(index, 0);
  constexpr int colors[] = {kBlue, kGreen, kRed, kOrange, kPink, kYellow-2};
  const static int ncolors = sizeof(colors) / sizeof(*colors);
  const int color = colors[iline % ncolors];
  for(auto& hit : hseed.hits()) {
    const double z = hit.pos().z();
    plotHitPhiZ(z, hseed.helix().circleAzimuth(z), color, iline);
  }
}

//-----------------------------------------------------------------------------
// Plot the phi-z hits associated with a sim particle
void AgnosticHelixFinderDiag::plotSimPhiZ(const SimParticle* sim, bool resolve, const int index) {
  if(!sim || !_data || !_data->chColl) return;

  const int sim_id = sim->id().asInt();
  if(_simInfo.find(sim_id) == _simInfo.end()) return;

  const int pdg = std::abs(sim->pdgId());
  if(pdg > 10000 || (pdg == 2212 && !_showProtons)) return;

  // Get the MC circle info
  float xC, yC, r, lambda, phi0;
  MCCircle(sim_id, xC, yC, r, lambda, phi0);
  if(_debugLevel > 2) printf(" Sim phi-z: r = %6.1f, lambda = %6.1f, phi 0 = %4.1f\n",
                             r, lambda, phi0);

  // Loop through the hits, drawing MC-matched ones
  constexpr int colors[] = {kBlue, kGreen, kRed, kOrange, kPink, kYellow-2};
  const static int ncolors = sizeof(colors) / sizeof(*colors);
  const int iline = std::max(index, 0);
  const int color = colors[iline % ncolors];

  const size_t nhits = _data->chColl->size();
  for(size_t index = 0; index < nhits; ++index) {
    const SimParticle* hit_sim = hitSim(index);
    if(hit_sim && hit_sim->id() == sim->id()) {
      const auto hit = _data->chColl->at(index);
      const double z = hit.pos().z();
      double phi = phi0 + z/lambda;
      if(!resolve) {
        phi = std::fmod(phi, 2.*M_PI);
        if(phi < 0.) phi += 2.*M_PI;
      }
      plotHitPhiZ(z, phi, color, iline, 2);
    }
  }
}

//-----------------------------------------------------------------------------
// Plot a fit circle in x-y
void AgnosticHelixFinderDiag::plotCircleXY() {
  if(!_data || !_data->tcHits || !_data->chColl) return;

  // check for MC info, plot below the reco data
  const auto sim = circleSimMatch();
  const int sim_id = (sim) ? sim->id().asInt() : -1;
  if(sim && _simInfo.find(sim_id) != _simInfo.end()) {
    float xC_mc, yC_mc, r_mc;
    MCCircle(sim_id, xC_mc, yC_mc, r_mc);
    const int pdg = std::abs(sim->pdgId());
    const int mc_color = MCColor(pdg);
    plotCircle(xC_mc, yC_mc, r_mc, mc_color, 2, "MCCircle"); // MC circle estimate

    // all of the MC hits
    for(size_t index = 0; index < _data->chColl->size(); ++index) {
      const auto hit_sim = hitSim(index);
      if(!hit_sim) continue;
      const int hit_sim_id = hit_sim->id().asInt();
      if(sim_id == hit_sim_id) {
        plotHitXY(_data->chColl->at(index), mc_color);
      }
    }
  }

  // Plot the reco circle
  ::LsqSums4 fitter = *(_data->circleFitter); // clone to not change the fitter
  const float xC = fitter.x0();
  const float yC = fitter.y0();
  const float r  = fitter.radius();
  plotCircle(xC, yC, r, kRed, 2, "RecoCircle"); // Reco circle

  // Plot each hit on the circle
  for(auto& tcHit : *(_data->tcHits)) {
    if(tcHit.inHelix || !tcHit.used) continue;
    const int index = tcHit.hitIndice;
    if(index < 0) continue;
    plotHitXY(_data->chColl->at(index), kGreen);
  }
}

//-----------------------------------------------------------------------------
// Plot a fit circle in x-y
void AgnosticHelixFinderDiag::plotCircleXY(const seedCircleInfo& info, int index) {
  if(!_data || !_data->tcHits || !_data->chColl) return;

  // Plot the reco circle
  const float xC = info.xC;
  const float yC = info.yC;
  const float r  = info.radius;
  if(_debugLevel > 2) printf("  %s: Plotting circle (x, y, r) = (%7.1f, %7.1f, %6.1f), index = %i\n",
                             __func__, xC, yC, r, index);
  plotCircle(xC, yC, r, kRed, 2, (index < 0) ? "RecoCircle" : Form("RecoCircle_%i", index)); // Reco circle

  // Plot each hit on the circle
  for(auto index : info.hits) {
    if(index < 0) continue;
    plotHitXY(_data->chColl->at(index), kGreen);
  }
}

//-----------------------------------------------------------------------------
// Plot the hits in phi-z for a reco circle
void AgnosticHelixFinderDiag::plotCirclePhiZ() {
  if(!_data || !_data->tcHits || !_data->chColl) return;

  // check for MC info, plot below the reco data
  const auto sim = circleSimMatch();
  if(sim) plotSimPhiZ(sim, false);

  // Plot the hits in phi-z
  ::LsqSums4 fitter = *(_data->circleFitter); // clone to not change the fitter
  const float xC = fitter.x0();
  const float yC = fitter.y0();

  // Plot each hit on the circle
  for(auto& tcHit : *(_data->tcHits)) {
    if(tcHit.inHelix || !tcHit.used) continue;
    const int index = tcHit.hitIndice;
    if(index < 0) continue;
    const auto& hit = _data->chColl->at(index);
    const float x = hit.pos().x() - xC;
    const float y = hit.pos().y() - yC;
    const float z = hit.pos().z();
    float phi = std::fmod(polyAtan2(y, x), 2.*M_PI);
    if(phi < 0.) phi += 2.*M_PI;
    plotHitPhiZ(z, phi, kGreen);
  }
}

//-----------------------------------------------------------------------------
// Plot a triplet in x-y
void AgnosticHelixFinderDiag::plotTripletXY(const tripletInfo& info, const int index, bool mc_triplet) {

  // Plot the circle from the triplet
  const float r  = info.radius;
  const float xC = info.xC;
  const float yC = info.yC;
  if(_debugLevel > 2) printf("    Triplet   %i: x = %6.1f, y = %6.1f, r = %6.1f\n", index, xC, yC, r);
  std::string title = (index >= 0) ? Form("TripletReco_%i", index) : "TripletReco";
  plotCircle(xC, yC, r, kRed, 2, title); // Reco triplet

  // Plot the triplet points used
  const int color = kGreen;
  plotHitXY(info.trip.i, color);
  plotHitXY(info.trip.j, color);
  plotHitXY(info.trip.k, color);

  // check for MC info
  const auto sim = tripletSimMatch(info.trip);
  const int sim_id = (sim) ? sim->id().asInt() : -1;
  if(sim && _simInfo.find(sim_id) != _simInfo.end()) {
    float xC_mc, yC_mc, r_mc;
    MCCircle(sim_id, xC_mc, yC_mc, r_mc);
    if(_debugLevel > 2) printf("    MC circle %i: x = %6.1f, y = %6.1f, r = %6.1f\n", index, xC_mc, yC_mc, r_mc);
    title = (index >= 0) ? Form("TripletMC_%i", index) : "TripletMC";
    plotCircle(xC_mc, yC_mc, r_mc, kBlue, 2, title); // MC helix estimate
  }
}

//-----------------------------------------------------------------------------
// Create axes for the XY plots
void AgnosticHelixFinderDiag::plotXYAxes() {
  // Use a graph to control the axes
  TGraph* axes = new TGraph();
  axes->AddPoint(0.,0.); // necessary for drawing
  axes->GetXaxis()->SetLimits(-900, 900);
  axes->GetYaxis()->SetRangeUser(-900, 900);
  axes->GetXaxis()->SetTitle("x (mm)");
  axes->GetYaxis()->SetTitle("y (mm)");
  axes->Draw("AP");
  _drawn_objects.push_back(axes);

  // Draw the tracker
  plotCircle(0., 0., 380., kBlack, 2, "tracker");
  plotCircle(0., 0., 700., kBlack, 2, "tracker", false);

  plotXYLegend(0);
}

//-----------------------------------------------------------------------------
// Create legend for the x-y plots
void AgnosticHelixFinderDiag::plotXYLegend(int) {
  TLegend* leg = new TLegend(0.16, 0.82, 0.80, 0.9);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  TGraph* g_r = new TGraph();
  g_r->SetLineColor(kRed);
  g_r->SetLineWidth(2);
  _drawn_objects.push_back(g_r);
  leg->AddEntry(g_r, "Reco", "L");
  TGraph* g_h = new TGraph();
  g_h->SetLineColor(kGreen);
  g_h->SetLineWidth(2);
  _drawn_objects.push_back(g_h);
  leg->AddEntry(g_h, "Used hits", "L");
  if(_simcol) {
    leg->SetNColumns(3);
    TGraph* g_m = new TGraph();
    g_m->SetLineColor(kBlue);
    g_m->SetLineWidth(2);
    leg->AddEntry(g_m, "MC", "L");
    _drawn_objects.push_back(g_m);
  } else {
    leg->SetNColumns(2);
  }
  leg->Draw();
  _drawn_objects.push_back(leg);
}

//-----------------------------------------------------------------------------
// Create axes for the phi-z plots
void AgnosticHelixFinderDiag::plotPhiZAxes(double phi_min, double phi_max) {
  // Use a graph to control the axes
  TGraph* axes = new TGraph();
  axes->AddPoint(0.,0.); // necessary for drawing
  axes->SetMarkerSize(0.);
  axes->GetXaxis()->SetLimits(-1800., 1800.);
  axes->GetYaxis()->SetRangeUser(phi_min, phi_max);
  axes->GetXaxis()->SetTitle("z (mm)");
  axes->GetYaxis()->SetTitle("#phi");
  axes->Draw("AP");
  _drawn_objects.push_back(axes);

  axes->GetXaxis()->SetLabelSize(0.08);
  axes->GetYaxis()->SetLabelSize(0.08);
  axes->GetXaxis()->SetTitleSize(0.10);
  axes->GetYaxis()->SetTitleSize(0.10);
  axes->GetYaxis()->SetTitleOffset(0.3);
  plotPhiZLegend(0);
}

//-----------------------------------------------------------------------------
// Create legend for the phi-z plots
void AgnosticHelixFinderDiag::plotPhiZLegend(int) {
  TLegend* leg = new TLegend(0.10, 0.8, 0.23, 0.9);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.10);
  TGraph* g_r = new TGraph();
  g_r->SetMarkerColor(kBlue);
  g_r->SetMarkerStyle(24);
  g_r->SetMarkerSize(2);
  _drawn_objects.push_back(g_r);
  leg->AddEntry(g_r, "Reco", "P");
  if(_simcol) {
    leg->SetNColumns(2);
    TGraph* g_m = new TGraph();
    g_m->SetMarkerColor(kBlue);
    g_m->SetMarkerStyle(2);
    g_m->SetMarkerSize(2);
    leg->AddEntry(g_m, "MC", "P");
    _drawn_objects.push_back(g_m);
  }
  leg->Draw();
  _drawn_objects.push_back(leg);
}

//-----------------------------------------------------------------------------
// Plot XY hits
void AgnosticHelixFinderDiag::plotHitsXY(const bool use_tc) {
  if(!_data || !_data->chColl) return;
  if(use_tc && !_data->tc) return;
  if(!use_tc && !_data->tcColl) return;

  const size_t ntc = (use_tc) ? 1 : _data->tcColl->size();
  for(size_t itc = 0; itc < ntc; ++itc) {
    const TimeCluster* tc = (use_tc) ? _data->tc : &_data->tcColl->at(itc);
    const size_t nhits = tc->hits().size();
    for(size_t ihit = 0; ihit < nhits; ++ihit) {
      const size_t index = tc->hits().at(ihit);
      const auto& hit = _data->chColl->at(index);
      if(hit.time() < _tMin) continue;
      const auto sim = hitSim(index);
      double pmc(0.); int pdg = 0; // for MC-based hit coloring, if MC is available
      if(sim) {
        pdg = std::abs(sim->pdgId());
        const int sim_id = sim->id().asInt();
        if(_simInfo.find(sim_id) != _simInfo.end()) {
          const auto& info = _simInfo[sim_id];
          pmc = std::max(info.hit_start_p_, info.hit_end_p_);
        }
        if(pdg == 2212 && !_showProtons) continue;
      }
      const int color = MCColor(pdg, pmc);
      plotHitXY(hit, color);
    }
  }
}

//-----------------------------------------------------------------------------
// Plot all relevant MC
void AgnosticHelixFinderDiag::plotMC(int stage, bool phiz) {
  // Relevant sim particles
  int index = 0;
  for(auto info_pair : _simInfo) {
    auto& info = info_pair.second;
    const double pmc = std::max(info.hit_start_p_, info.hit_end_p_);
    const double tmc = std::max(info.hit_start_t_, info.hit_end_t_);
    if(info.nhits_ > 10 && pmc > _pMin && tmc > _tMin) {
      const auto sim = simByID(_simcol, info.id_);
      const int pdg = std::abs(sim->pdgId());
      if((pdg < 10000) && (pdg != 2212 || _showProtons)) {
        if(stage == kHelix) { // helix canvas
          if(phiz) {
            plotSimPhiZ(sim, true, index);
          } else {
            float xC_mc, yC_mc, r_mc;
            MCCircle(info.id_, xC_mc, yC_mc, r_mc);
            if(_debugLevel > 2) printf(" Relevant MC: pdg = %4i, pmc = %6.1f, x(C) = %7.1f, y(C) = %7.1f, r = %6.1f\n",
                                     pdg, pmc, xC_mc, yC_mc, r_mc);
            const int mc_color = MCColor(pdg, pmc);
            plotCircle(xC_mc, yC_mc, r_mc, mc_color, 2, "MC_Sim"); // MC helix estimate
          }
        } else if(stage == kTriplet && !phiz) { // triplets
          if(_simTriplets.find(info.id_) != _simTriplets.end()) {
            plotTripletXY(_simTriplets[info.id_], index, true);
          }
        } else if(stage == kCircle && !phiz) { // seed circles
          if(_simCircles.find(info.id_) != _simCircles.end()) {
            plotCircleXY(_simCircles[info.id_], index);
          } else {
            if(_debugLevel > 0) printf("--> Missing MC circle for ID %i, index %i\n", info.id_, index);
          }
        }
        ++index;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Plot the hits in an event
void AgnosticHelixFinderDiag::plotBeginStage(const bool use_tc) {
  if(_debugLevel > 0) printf("  %12s: Plotting the event, use timecluster = %o\n",
                             __func__, use_tc);
  auto c = _c_hlx;
  if(!beginPlot(c)) return;
  std::string title("Hits");
  plotXYAxes();
  plotHitsXY(use_tc);
  plotMC(kHelix);
  endPlot(c, title);
}

//-----------------------------------------------------------------------------
// Plot the triplet(s)
void AgnosticHelixFinderDiag::plotTripletStage(const bool mc_triplets) {
  if(_debugLevel > 0) printf("  %12s: Plotting the triplet circle, MC = %o\n",
                             __func__, mc_triplets);
  auto c = _c_trp;
  if(!beginPlot(c)) return;
  std::string title("Triplet");
  plotXYAxes();
  if(mc_triplets) {  // all stored triplets best matched to relevant sims
    plotMC(kTriplet);
    title = "Triplets best-matched to MC";
  }
  else plotTripletXY(_data->tripInfo); // just the current triplet
  endPlot(c, title);
}

//-----------------------------------------------------------------------------
// Plot the seed circle(s)
void AgnosticHelixFinderDiag::plotCircleStage(const bool mc_circles) {
  if(_debugLevel > 0) printf("  %12s: Plotting the seed circle, MC = %o\n",
                             __func__, mc_circles);
  auto c = (mc_circles) ? _c_trp : _c_hlx;
  if(!beginPlot(c)) return;
  std::string title("Seed Circle");
  plotXYAxes();
  if(mc_circles) {
    plotMC(kCircle);
    title = "Seed circles matched to MC";
  } else {
    plotCircleXY();
  }

  endPlot(c, title);
}

//-----------------------------------------------------------------------------
// Plot the seed circle(s)
void AgnosticHelixFinderDiag::plotSegmentStage(bool resolve, bool seed_circle) {
  if(_debugLevel > 0) printf("  %12s: Plotting the segments, resolve = %o\n",
                             __func__, resolve);
  auto c = _c_seg;
  if(!beginPlot(c)) return;
  std::string title("Line segment #phi-z");

  // Draw the axes
  if(seed_circle || !resolve) plotPhiZAxes(-0.2, 6.5);
  else { // determine the phi range needed
    double phi_min(-2.*M_PI), phi_max(2.*M_PI); // default to -2pi - 2pi
    for(size_t iline = 0; iline < _data->seedPhiLines->size(); ++iline) {
      const auto& line = _data->seedPhiLines->at(iline);
      const double zmin  = line.zMin;
      const double zmax  = line.zMax;
      // clone the fitter to avoid changing its state
      ::LsqSums2 fitter = line.fitter;
      const double slope = fitter.dydx();
      const double phi0  = fitter.y0();;
      const double phi1  = slope * zmin + phi0;
      const double phi2  = slope * zmax + phi0;
      phi_min = std::min(phi_min, std::min(phi1, phi2));
      phi_max = std::max(phi_max, std::max(phi1, phi2));
    }
    const double buffer = 0.1*(phi_max - phi_min);
    plotPhiZAxes(phi_min - buffer, phi_max + 2.*buffer);
  }

  // Draw the info
  if(seed_circle) {
    title = "Seed Circle #phi-z";
    plotCirclePhiZ();
  } else {
    for(size_t iline = 0; iline < _data->seedPhiLines->size(); ++iline) {
      plotLinePhiZ(_data->seedPhiLines->at(iline), false, iline);
    }
  }

  endPlot(c, title);
}

//-----------------------------------------------------------------------------
// Plot the helix (helices) in XY
void AgnosticHelixFinderDiag::plotHelixStageXY(int stage) {
  if(_debugLevel > 0) printf("  %12s: Plotting the helix XY circle, stage = %i\n",
                             __func__, stage);
  auto c = _c_hlx;
  if(!beginPlot(c)) return;
  std::string title("Helices");
  plotXYAxes();
  if(stage == kHelix) { // Current helix
    title = "Helix";
    plotHelixXY(*(_data->hseed));
  } else { // All helices
    const bool use_tc = stage == kTimeCluster;
    if(use_tc) title = "Helices (time cluster)";
    plotHitsXY(use_tc); // draw (relevant) hits
    plotMC(kHelix);
    for(size_t index = 0; index < _data->helixSeedData.size(); ++index) {
      plotHelixXY(_data->helixSeedData.at(index).seed, index);
    }
  }

  endPlot(c, title);
}

//-----------------------------------------------------------------------------
// Plot the helix (helices) in phi-z
void AgnosticHelixFinderDiag::plotHelixStagePhiZ(int stage) {
  if(_debugLevel > 0) printf("  %12s: Plotting the helix phi-z line, stage = %i\n",
                             __func__, stage);
  auto c = _c_seg;
  if(!beginPlot(c)) return;
  std::string title("Helices #phi-z");

  // determin the phi range needed
  double phi_min(-2.*M_PI), phi_max(2.*M_PI);
  if(stage == kHelix) {
    const auto& seed = *(_data->hseed);
    const double phi_1 = seed.helix().circleAzimuth(-1600.);
    const double phi_2 = seed.helix().circleAzimuth( 1600.);
    phi_min = std::min(phi_min, std::min(phi_1, phi_2));
    phi_max = std::max(phi_max, std::max(phi_1, phi_2));
  } else {
    for(size_t index = 0; index < _data->helixSeedData.size(); ++index) {
      const auto& seed = _data->helixSeedData.at(index).seed;
      const double phi_1 = seed.helix().circleAzimuth(-1600.);
      const double phi_2 = seed.helix().circleAzimuth( 1600.);
      phi_min = std::min(phi_min, std::min(phi_1, phi_2));
      phi_max = std::max(phi_max, std::max(phi_1, phi_2));
    }
  }
  const double buffer = 0.1*(phi_max - phi_min);
  plotPhiZAxes(phi_min - buffer, phi_max + 2.*buffer);

  // Plot the info
  if(stage == kHelix) { // Current helix
    title = "Helix #phi-z";
    const auto& seed = *(_data->hseed);
    const auto sim = helixSimMatch(seed);
    if(sim) plotSimPhiZ(sim, true);
    plotHelixPhiZ(seed);
  } else { // All helices
    plotMC(kHelix, true);
    for(size_t index = 0; index < _data->helixSeedData.size(); ++index) {
      plotHelixPhiZ(_data->helixSeedData.at(index).seed, index);
    }
  }

  endPlot(c, title);
}

//-----------------------------------------------------------------------------
// Plot the 3D event: Not fully implemented
void AgnosticHelixFinderDiag::plotTotal(int option) {
  if(_debugLevel > -1) printf("  %12s: Plotting with option %i\n", __func__, option);
  if(!_data || !_data->tcHits || !_data->chColl) return;

  // Retrieve the relevant canvas and switch to this pad
  auto c = _c_3d;
  if(!c) {
    printf("AgnosticHelixFinderDiag::%s: Undefined option %i\n", __func__, option);
    return;
  }
  c->cd();
  c->Draw();

  // Use a graph to control the axes
  TH3F* frame3d = new TH3F("frame3d","frame3d",
                           10, -1600., 1600.,
                           10, -900, 900,
                           10, -900, 900);
  _drawn_objects.push_back(frame3d);

  // Make a graph of all hits
  auto g_hits = new TGraph2D();
  for(const auto& hit : *(_data->chColl)) {
    const auto pos = hit.pos();
    g_hits->AddPoint(pos.z(), pos.x(), pos.y());
  }
  g_hits->SetMarkerStyle(20);
  g_hits->SetMarkerColor(kRed);
  _drawn_objects.push_back(g_hits);

  // Draw the "tracker" -- just use points for simplicity
  TGraph2D* g_tracker = new TGraph2D();
  const int nphi = 40; const int nz = 4;
  const double rmin(380.), rmax(700.), zmin(-1500.), zmax(1500.);
  for(int iphi = 0; iphi <= nphi; ++iphi) {
    for(int iz = 0; iz <= nz; ++iz) {
      const double phi = 2.*M_PI*iphi/nphi;
      const double z = zmin + (zmax - zmin)*iz/nz;
      g_tracker->AddPoint(z, rmin*std::cos(phi), rmin*std::sin(phi));
      g_tracker->AddPoint(z, rmax*std::cos(phi), rmax*std::sin(phi));
    }
  }
  g_tracker->SetMarkerStyle(6);
  g_tracker->SetMarkerColor(kBlack);
  _drawn_objects.push_back(g_tracker);

  // FIXME: Currently 3D drawing returns warnings due to ROOT issue and art upgrades these to errors
  if(false) {
    frame3d->Draw();
    g_tracker->Draw("P SAME");
    g_hits->Draw("P SAME");
  } else { // only draw the hits for now
    g_hits->Draw("P SAME");
  }

  plotLabel(c->GetLeftMargin(), c->GetTopMargin());
  c->Modified(); c->Update();
}
