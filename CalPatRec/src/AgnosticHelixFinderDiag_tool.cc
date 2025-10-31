// Tool to analzye/debug Agnostic Helix Finder processing
// Original author: Matthew Stortini

#include "Offline/CalPatRec/inc/AgnosticHelixFinderDiag.hh"

using namespace mu2e;

//-----------------------------------------------------------------------------
AgnosticHelixFinderDiag::AgnosticHelixFinderDiag(const fhicl::Table<mu2e::AgnosticHelixFinderTypes::Config>& config) :
  _simTag(config().simTag())
  , _digTag(config().digTag())
  , _display(config().display())
  , _display3D(config().display3D())
  , _showProtons(config().showProtons())
  , _pMin(config().pMin())
  , _tMin(config().tMin())
  , _debugLevel(config().debug())
{
  if(gROOT->IsBatch()) _display = false;
}

//-----------------------------------------------------------------------------
AgnosticHelixFinderDiag::~AgnosticHelixFinderDiag() {}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir) {

  Hist->nHelices = Dir->make<TH1F>("nHelices", "number of helices found per event", 30, 0.0, 30.0);
  Hist->nTimeClusters =
    Dir->make<TH1F>("nTimeClusters", "number of time clusters per event", 100, 0.0, 100.0);
  Hist->MC_nSim = Dir->make<TH1F>("MC_nSim", "N(relevant sims)", 30, 0.0, 30.0);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir) {

  Hist->nHelicesPerTC =
    Dir->make<TH1F>("nHelicesPerTC", "number of helices found per TC", 30, 0, 30.0);
  Hist->nComboHitsPerTC =
    Dir->make<TH1F>("nComboHitsPerTC", "number of combo hits per TC", 200, 0, 200.0);
  Hist->nStrawHitsPerTC =
    Dir->make<TH1F>("nStrawHitsPerTC", "number of straw hits per TC", 300, 0, 300.0);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookHelixSeedHistograms(HelixSeedHists* Hist, art::TFileDirectory* Dir) {

  Hist->eDepAvg        = Dir->make<TH1F>("eDepAvg"         , "average combo hit E_{dep}" , 1000,    0.,  0.01);
  Hist->p              = Dir->make<TH1F>("p"               , "Reco momentum"             ,  400,    0.,  400.);
  Hist->t0             = Dir->make<TH1F>("t0"              , "Reco t_{0}"                ,  200,    0., 2000.);
  Hist->radius         = Dir->make<TH1F>("radius"          , "Reco radius"               ,  200,    0., 1000.);
  Hist->lambda         = Dir->make<TH1F>("lambda"          , "Reco #lambda"              ,  200,    0., 2000.);
  Hist->d0             = Dir->make<TH1F>("d0"              , "Reco d_{0}"                ,  200, -600.,  600.);
  Hist->hits           = Dir->make<TH1F>("hits"            , "N(hits)"                   ,  200,    0.,   200.);

  Hist->MC_PDG         = Dir->make<TH1F>("MC_PDG"          , "MC PDG ID"                 , 5000, -2500, 2500.); // FIXME: Store PDG name, not ID
  Hist->MC_p           = Dir->make<TH1F>("MC_p"            , "MC momentum"               ,  400,    0.,  400.);
  Hist->MC_dp          = Dir->make<TH1F>("MC_dp"           , "p_{reco} - p_{MC}"         ,  400,  -50.,   50.);
  Hist->MC_purity      = Dir->make<TH1F>("MC_purity"       , "N(correct hits)/N(hits)"   ,  110,    0.,   1.1);
  Hist->MC_hitFrac     = Dir->make<TH1F>("MC_hitFrac"      , "N(MC hits used)/N(MC hits)",  110,    0.,   1.1);
  Hist->MC_hits        = Dir->make<TH1F>("MC_hits"         , "N(MC hits)"                ,  200,    0.,   200.);
  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookLineInfoHistograms(LineInfoHists* Hist, art::TFileDirectory* Dir) {

  Hist->nHitsRatioPerHS =
    Dir->make<TH1F>("nHitsRatioPerHS", "ratio of circle fit hits to phi fit hits", 1000, 0, 5.0);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookLineSegmentHistograms(LineSegmentHists* Hist, art::TFileDirectory* Dir) {

  Hist->chi2dof = Dir->make<TH1F>("chi2dof", "chi2dof of line segments", 300, 0, 30.0);
  Hist->maxHitGap = Dir->make<TH1F>("maxHitGap", "max dZ between adjacent points", 600, 0, 600.0);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookTripletHistograms(TripletHists* Hist, art::TFileDirectory* Dir) {

  Hist->radius         = Dir->make<TH1F>("radius"          , "Reco radius"               ,  200,    0., 1000.);
  Hist->rmax           = Dir->make<TH1F>("rmax"            , "Reco R(max)"               ,  300,    0., 2000.);
  Hist->d0             = Dir->make<TH1F>("d0"              , "Reco d_{0}"                ,  300,-1000., 1000.);
  Hist->dz12           = Dir->make<TH1F>("dz12"            , "#Deltaz_{12}"              ,  200,    0., 3000.);
  Hist->dz23           = Dir->make<TH1F>("dz23"            , "#Deltaz_{23}"              ,  200,    0., 3000.);
  Hist->dz13           = Dir->make<TH1F>("dz13"            , "#Deltaz_{13}"              ,  200,    0., 3000.);
  Hist->d12            = Dir->make<TH1F>("d12"             , "(x,y) distance_{12}"       ,  200,    0., 2000.);
  Hist->d23            = Dir->make<TH1F>("d23"             , "(x,y) distance_{23}"       ,  200,    0., 2000.);
  Hist->d13            = Dir->make<TH1F>("d13"             , "(x,y) distance_{13}"       ,  200,    0., 2000.);
  Hist->dphi12         = Dir->make<TH1F>("dphi12"          , "#Delta#phi_{12}"           ,  100,  -3.2,   3.2);
  Hist->dphi23         = Dir->make<TH1F>("dphi23"          , "#Delta#phi_{23}"           ,  100,  -3.2,   3.2);
  Hist->dphi13         = Dir->make<TH1F>("dphi13"          , "#Delta#phi_{13}"           ,  100,  -3.2,   3.2);

  Hist->MC_PDG_1       = Dir->make<TH1F>("MC_PDG_1"        , "MC PDG ID(1)"              , 5000, -2500,  2500);
  Hist->MC_PDG_2       = Dir->make<TH1F>("MC_PDG_2"        , "MC PDG ID(2)"              , 5000, -2500,  2500);
  Hist->MC_PDG_3       = Dir->make<TH1F>("MC_PDG_3"        , "MC PDG ID(3)"              , 5000, -2500,  2500);
  Hist->MC_p           = Dir->make<TH1F>("MC_p"            , "MC momentum"               ,  400,     0,  400);
  Hist->MC_pt          = Dir->make<TH1F>("MC_pt"           , "MC p_{T}"                  ,  400,     0,  400);
  Hist->MC_dr          = Dir->make<TH1F>("MC_dr"           , "R(reco) - R(MC)"           ,  200,  -500,  500);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
  char folder_name[20];
  TH1::AddDirectory(0);

  //-----------------------------------------------------------------------------
  // book event histograms
  //-----------------------------------------------------------------------------
  std::string ev_titles[kNEventHistsSets];
  ev_titles[0] = "All events";
  for (int i = 0; i < kNEventHistsSets; i++) {
    if(ev_titles[i] == "") continue;
    sprintf(folder_name, "evt_%i", i);
    art::TFileDirectory tfdir = Tfs->mkdir(folder_name, ev_titles[i].c_str());
    _hist._eventHists[i] = new EventHists;
    bookEventHistograms(_hist._eventHists[i], &tfdir);
  }

  //-----------------------------------------------------------------------------
  // book time cluster histograms
  //-----------------------------------------------------------------------------
  std::string tc_titles[kNTimeClusterHistsSets];
  tc_titles[0] = "All time clusters";
  for (int i = 0; i < kNTimeClusterHistsSets; i++) {
    if(tc_titles[i] == "") continue;
    sprintf(folder_name, "tcl_%i", i);
    art::TFileDirectory tfdir = Tfs->mkdir(folder_name, tc_titles[i].c_str());
    _hist._timeClusterHists[i] = new TimeClusterHists;
    bookTimeClusterHistograms(_hist._timeClusterHists[i], &tfdir);
  }

  //-----------------------------------------------------------------------------
  // book helix seed histograms
  //-----------------------------------------------------------------------------
  std::string hs_titles[kNHelixSeedHistsSets];
  hs_titles[0] = "All helices";
  hs_titles[1] = "Accepted helices";
  for (int i = 0; i < kNHelixSeedHistsSets; i++) {
    if(hs_titles[i] == "") continue;
    sprintf(folder_name, "hs_%i", i);
    art::TFileDirectory tfdir = Tfs->mkdir(folder_name, hs_titles[i].c_str());
    _hist._helixSeedHists[i] = new HelixSeedHists;
    bookHelixSeedHistograms(_hist._helixSeedHists[i], &tfdir);
  }

  //-----------------------------------------------------------------------------
  // book lineInfo histograms
  //-----------------------------------------------------------------------------
  std::string ln_titles[kNLineInfoHistsSets];
  ln_titles[0] = "All lines";
  for (int i = 0; i < kNLineInfoHistsSets; i++) {
    if(ln_titles[i] == "") continue;
    sprintf(folder_name, "li_%i", i);
    art::TFileDirectory tfdir = Tfs->mkdir(folder_name, ln_titles[i].c_str());
    _hist._lineInfoHists[i] = new LineInfoHists;
    bookLineInfoHistograms(_hist._lineInfoHists[i], &tfdir);
  }

  //-----------------------------------------------------------------------------
  // book line segment histograms
  //-----------------------------------------------------------------------------
  std::string ls_titles[kNLineSegmentHistsSets];
  ls_titles[0] = "All line segments";
  for (int i = 0; i < kNLineSegmentHistsSets; i++) {
    if(ls_titles[i] == "") continue;
    sprintf(folder_name, "ls_%i", i);
    art::TFileDirectory tfdir = Tfs->mkdir(folder_name, ls_titles[i].c_str());
    _hist._lineSegmentHists[i] = new LineSegmentHists;
    bookLineSegmentHistograms(_hist._lineSegmentHists[i], &tfdir);
  }

  //-----------------------------------------------------------------------------
  // book triplet histograms
  //-----------------------------------------------------------------------------
  std::string tp_titles[kNTripletHistsSets];
  tp_titles[0] = "All triplets";
  tp_titles[1] = "Triplets with lines";
  tp_titles[2] = "Triplets with helices";
  tp_titles[3] = "Accepted triplets";
  for (int i = 0; i < kNTripletHistsSets; i++) {
    if(tp_titles[i] == "") continue;
    sprintf(folder_name, "tp_%i", i);
    art::TFileDirectory tfdir = Tfs->mkdir(folder_name, tp_titles[i].c_str());
    _hist._tripletHists[i] = new TripletHists;
    bookTripletHistograms(_hist._tripletHists[i], &tfdir);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillEventHistograms(EventHists* Hist, diagInfo* Data) {

  Hist->nHelices->Fill(Data->nHelices);
  Hist->nTimeClusters->Fill(Data->nTimeClusters);

  // MC info
  int nsim(0);
  for(auto& sim_pair : _simInfo) {
    auto& sim = sim_pair.second;
    if(sim.nhits_ > 10 && (sim.hit_start_p_ > 20. || sim.hit_end_p_ > 20.)) ++nsim;
  }
  Hist->MC_nSim->Fill(nsim);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillTimeClusterHistograms(TimeClusterHists* Hist, diagInfo* Data,
                                                       int loopIndex) {

  // fill per tc info
  Hist->nHelicesPerTC->Fill(Data->timeClusterData.at(loopIndex).nHelices);
  Hist->nComboHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nComboHits);
  Hist->nStrawHitsPerTC->Fill(Data->timeClusterData.at(loopIndex).nStrawHits);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillHelixSeedHistograms(HelixSeedHists* Hist, const hsInfo& info) {

  // Retrieve the info
  const HelixSeed& seed = info.seed;
  if(seed.hits().empty()) return 1; // empty helix seed

  const float bz0 = info.bz0;
  constexpr float mmMeV = 3./10.; // convert to MeV/c
  const float edep  = seed.hits().eDepAvg();
  const float p     = mmMeV * bz0 * seed.helix().momentum();
  const float t     = seed.t0().t0();
  const float r     = seed.helix().radius();
  const float l     = std::fabs(seed.helix().lambda());
  const float d0    = seed.helix().rcent() - r; // no helicity sign, for simplicity
  const int   nhits = seed.hits().nStrawHits();

  // Reco info
  Hist->eDepAvg->Fill(edep);
  Hist->p      ->Fill(p);
  Hist->t0     ->Fill(t);
  Hist->radius ->Fill(r);
  Hist->lambda ->Fill(l);
  Hist->d0     ->Fill(d0);
  Hist->hits   ->Fill(nhits);

  // MC info, if available
  double pmc; int mc_hits;
  const auto sim = helixSimMatch(seed, pmc, mc_hits);
  const int pdg = (sim) ? sim->pdgId() : 0;

  Hist->MC_PDG ->Fill(pdg);
  Hist->MC_p   ->Fill(pmc);
  Hist->MC_dp  ->Fill(p - pmc);
  Hist->MC_purity  ->Fill(MCHitPurity(seed));
  Hist->MC_hitFrac ->Fill(MCHitFraction(seed));
  Hist->MC_hits    ->Fill(mc_hits);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillLineInfoHistograms(LineInfoHists* Hist, diagInfo* Data,
                                                    int loopIndex) {

  // fill per line info
  Hist->nHitsRatioPerHS->Fill(Data->lineInfoData.at(loopIndex).nHitsRatio);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillLineSegmentHistograms(LineSegmentHists* Hist, diagInfo* Data,
                                                       int loopIndex) {

  // fill per line segment info
  Hist->chi2dof->Fill(Data->lineSegmentData.at(loopIndex).chi2dof);
  Hist->maxHitGap->Fill(Data->lineSegmentData.at(loopIndex).maxHitGap);

  return 0;
}

//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillTripletHistograms(TripletHists* Hist, const tripletInfo& info) {

  const auto& trip = info.trip;
  const float r  = info.radius;
  const float xC = info.xC;
  const float yC = info.yC;
  const float rC = std::sqrt(xC*xC + yC*yC);
  const float rmax = r + rC;
  const float d0 = rC - r;
  const float dz12 = std::fabs(trip.i.pos.z() - trip.j.pos.z());
  const float dz23 = std::fabs(trip.j.pos.z() - trip.k.pos.z());
  const float dz13 = std::fabs(trip.i.pos.z() - trip.k.pos.z());
  const float d12  = std::sqrt((trip.i.pos - trip.j.pos).Perp2());
  const float d23  = std::sqrt((trip.j.pos - trip.k.pos).Perp2());
  const float d13  = std::sqrt((trip.i.pos - trip.k.pos).Perp2());

  const float phi1 = polyAtan2(trip.i.pos.x() - xC, trip.i.pos.y() - yC);
  const float phi2 = polyAtan2(trip.j.pos.x() - xC, trip.j.pos.y() - yC);
  const float phi3 = polyAtan2(trip.k.pos.x() - xC, trip.k.pos.y() - yC);
  float dphi12 = phi1 - phi2;
  float dphi23 = phi2 - phi3;
  float dphi13 = phi1 - phi3;
  if(dphi12 < -M_PI) dphi12 += 2.*M_PI;
  if(dphi12 >  M_PI) dphi12 -= 2.*M_PI;
  if(dphi23 < -M_PI) dphi23 += 2.*M_PI;
  if(dphi23 >  M_PI) dphi23 -= 2.*M_PI;
  if(dphi13 < -M_PI) dphi13 += 2.*M_PI;
  if(dphi13 >  M_PI) dphi13 -= 2.*M_PI;

  Hist->radius         ->Fill(r);
  Hist->rmax           ->Fill(rmax);
  Hist->d0             ->Fill(d0);
  Hist->dz12           ->Fill(dz12);
  Hist->dz23           ->Fill(dz23);
  Hist->dz13           ->Fill(dz13);
  Hist->d12            ->Fill(d12);
  Hist->d23            ->Fill(d23);
  Hist->d13            ->Fill(d13);
  Hist->dphi12         ->Fill(dphi12);
  Hist->dphi23         ->Fill(dphi23);
  Hist->dphi13         ->Fill(dphi13);

  // MC info
  const SimParticle* sim_1 = hitSim(trip.i.hitIndice);
  const SimParticle* sim_2 = hitSim(trip.j.hitIndice);
  const SimParticle* sim_3 = hitSim(trip.k.hitIndice);

  const SimParticle* best_sim = tripletSimMatch(trip);
  const int best_sim_id = (best_sim) ? best_sim->id().asInt() : -1;
  const auto* sim_info = (best_sim_id >= 0 && _simInfo.find(best_sim_id) != _simInfo.end()) ? &_simInfo[best_sim_id] : nullptr;

  const float p_mc  = (sim_info) ? 0.5*(sim_info->hit_start_p_  + sim_info->hit_end_p_ ) : 0.;
  const float pt_mc = (sim_info) ? 0.5*(sim_info->hit_start_pt_ + sim_info->hit_end_pt_) : 0.;

  float xC_mc, yC_mc, r_mc;
  MCCircle(best_sim_id, xC_mc, yC_mc, r_mc);

  Hist->MC_PDG_1->Fill((sim_1) ? sim_1->pdgId() : 0);
  Hist->MC_PDG_2->Fill((sim_2) ? sim_2->pdgId() : 0);
  Hist->MC_PDG_3->Fill((sim_3) ? sim_3->pdgId() : 0);
  Hist->MC_p->Fill(p_mc);
  Hist->MC_pt->Fill(pt_mc);
  Hist->MC_dr->Fill(r - r_mc);

  return 0;
}

//-----------------------------------------------------------------------------
// Mode used to differentiate steps in the helix search algorithm
//-----------------------------------------------------------------------------
int AgnosticHelixFinderDiag::fillHistograms(void* Data, int Mode) {

  if(!Data) return -1;
  _data = static_cast<diagInfo*>(Data);

  //-----------------------------------------------------------------------------
  // Initialize event info
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kBegin) {
    if(_debugLevel > 1) printf("[AgnosticHelixFinderDiag::%s: Processing the BEGIN stage\n", __func__);
    if(_data->event) { // the art::Event object is available
      // retrieve the collections
      art::Handle<SimParticleCollection> simH;
      art::Handle<StrawDigiMCCollection> digH;
      if(_simTag != "") _data->event->getByLabel<SimParticleCollection>(_simTag, simH);
      if(_digTag != "") _data->event->getByLabel<StrawDigiMCCollection>(_digTag, digH);
      _simcol = (simH.isValid()) ? simH.product() : nullptr;
      _digcol = (digH.isValid()) ? digH.product() : nullptr;
    } else {
      _simcol = nullptr;
      _digcol = nullptr;
    }
    initSimInfo(_simcol, _digcol);
    if(_display) {
      initDisplay();
      if(_data->diagLevel > 4) {
        plotBeginStage();
        runDisplay();
      }
    }
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Triplet stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kTriplet) {
    if(_display && _data->diagLevel > 5) {
      std::cout << "Triplet stage: Loop condition " << ConditionName(_data->loopCondition)
                << std::endl;
      plotTripletStage();
      plotBeginStage(true);
      runDisplay();
    }
    fillTripletHistograms(_hist._tripletHists[0], _data->tripInfo);
    checkSimTriplets();
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Circle stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kCircle) {
    if(_display && _data->diagLevel > 4) {
      std::cout << "Circle stage: Loop condition " << ConditionName(_data->loopCondition)
                << std::endl;
      plotTripletStage();
      plotCircleStage();
      plotSegmentStage(false, true); // circle in phi-z
      runDisplay();
    }
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Line segment stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kSegments) {
    if(_display && _data->diagLevel > 4) {
      std::cout << "Segment stage: Loop condition " << ConditionName(_data->loopCondition)
                << std::endl;
      plotTripletStage();
      plotSegmentStage(false);
      runDisplay();
    }
    fillTripletHistograms(_hist._tripletHists[1], _data->tripInfo);
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Line stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kLine) {
    if(_display && _data->diagLevel > 4) {
      std::cout << "Line stage: Loop condition " << ConditionName(_data->loopCondition)
                << std::endl;
      plotTripletStage();
      plotSegmentStage(true);
      runDisplay();
    }
    fillTripletHistograms(_hist._tripletHists[1], _data->tripInfo);
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Helix stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kHelix) {
    if(_data->hseed) {
      hsInfo info(_data->bz0, *(_data->hseed));
      if(_display && _data->diagLevel > 3) {
        std::cout << "Helix stage: Loop condition " << ConditionName(_data->loopCondition)
                  << std::endl;
        plotTripletStage();
        plotHelixStageXY(kHelix);
        plotHelixStagePhiZ(kHelix);
        runDisplay();
      }
      fillHelixSeedHistograms(_hist._helixSeedHists[0], info);
      fillTripletHistograms(_hist._tripletHists[2], _data->tripInfo);
    }
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Accepted Helix stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kFinal) {
    if(_data->hseed) {
      // only display at this level, where per-helix wasn't displayed but accepted helices are each plotted
      if(_display && _data->diagLevel == 3) {
        std::cout << "Helix stage: Loop condition " << ConditionName(_data->loopCondition)
                  << std::endl;
        plotTripletStage();
        plotHelixStageXY(kHelix);
        plotHelixStagePhiZ(kHelix);
        runDisplay();
      }
      hsInfo info(_data->bz0, *(_data->hseed));
      fillHelixSeedHistograms(_hist._helixSeedHists[1], info);
      fillTripletHistograms(_hist._tripletHists[3], _data->tripInfo);
    }
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Time cluster stage histograms
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kTimeCluster) {
    if(_display && _data->diagLevel > 1) {
      plotTripletStage(true);
      plotHelixStageXY(kTimeCluster);
      plotHelixStagePhiZ(kTimeCluster);
      runDisplay();
    }
    return 0;
  }

  //-----------------------------------------------------------------------------
  // Final histograms, only filled at event end
  //-----------------------------------------------------------------------------

  if(Mode == DIAG::kEnd) {
    if(_display) {
      plotTripletStage(true);
      plotHelixStageXY(kEnd);
      plotHelixStagePhiZ(kEnd);
      if(_display3D) plotTotal(0);
      runDisplay();
    }

    //-----------------------------------------------------------------------------
    // fill event histograms
    //-----------------------------------------------------------------------------
    fillEventHistograms(_hist._eventHists[0], _data);

    //-----------------------------------------------------------------------------
    // fill time cluster histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->timeClusterData.size(); i++) {
      fillTimeClusterHistograms(_hist._timeClusterHists[0], _data, i);
    }

    //-----------------------------------------------------------------------------
    // fill lineInfo seed histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->lineInfoData.size(); i++) {
      fillLineInfoHistograms(_hist._lineInfoHists[0], _data, i);
    }


    //-----------------------------------------------------------------------------
    // fill line segment histograms
    //-----------------------------------------------------------------------------
    for (int i = 0; i < (int)_data->lineSegmentData.size(); i++) {
      fillLineSegmentHistograms(_hist._lineSegmentHists[0], _data, i);
    }
  }

  return 0;
}

//--------------------------------------------------------------------------------------
void AgnosticHelixFinderDiag::checkSimTriplets() {
  if(!_data || !_simcol) return;

  // See if the current triplet is the best match the sim hits in it
  const auto& tripInfo = _data->tripInfo;

  // Check if this triplet is better for each sim
  for(auto& info_pair : _simInfo) {
    const auto& info = info_pair.second;
    if(info.nhits_ == 0) continue;
    const unsigned sim_id = info.id_;

    // Check if no triplet has been found
    if(_simTriplets.find(sim_id) == _simTriplets.end()) {
      _simTriplets[sim_id] = tripInfo;
      continue;
    }

    // If one has been found, determine which is better
    float x, y, r;
    MCCircle(sim_id, x, y, r);

    const float dr = tripInfo.radius - r;
    const float dx = tripInfo.xC - x;
    const float dy = tripInfo.yC - y;
    const float d = dr*dr + dx*dx + dy*dy;

    const auto& curr = _simTriplets[sim_id];
    const float dr_curr = curr.radius - r;
    const float dx_curr = curr.xC - x;
    const float dy_curr = curr.yC - y;
    const float d_curr = dr_curr*dr_curr + dx_curr*dx_curr + dy_curr*dy_curr;

    if(d < d_curr) _simTriplets[sim_id] = tripInfo;
  }
}

//--------------------------------------------------------------------------------------
float AgnosticHelixFinderDiag::MCHitPurity(const HelixSeed& hlx) {
  if(hlx.hits().empty()) return -1.f;
  int mc_hits;
  const auto sim = helixSimMatch(hlx, mc_hits);
  if(!sim) return -1.f;
  if(mc_hits < 0) return -1.f;
  const size_t nhits = hlx.hits().nStrawHits();
  return float(mc_hits) / nhits;
}

//--------------------------------------------------------------------------------------
float AgnosticHelixFinderDiag::MCHitFraction(const HelixSeed& hlx) {
  if(hlx.hits().empty()) return -1.f;
  int mc_hits_hlx;
  const auto sim = helixSimMatch(hlx, mc_hits_hlx);
  if(!sim) return -1.f;
  if(mc_hits_hlx < 0) return -1.f;
  const int mc_hits = MCHits(sim);
  if(mc_hits <= 0) return -1.f;
  return float(mc_hits_hlx) / mc_hits;
}

//--------------------------------------------------------------------------------------
void AgnosticHelixFinderDiag::initSimInfo(const SimParticleCollection* simcol,
                                          const StrawDigiMCCollection* digcol) {
  _simInfo.clear(); // clear the info from the previous event
  _simTriplets.clear();
  if(!simcol || !digcol) return;

  // Count the number of hits for each sim particle
  std::map<unsigned int,unsigned> pmap;
  std::map<unsigned int, Sim_t> info_map;
  if(_debugLevel > 5) printf("AgnosticHelixFinderDig::%s:\n", __func__);
  for(auto const& digi : *digcol) {
    // do not inspect truth info for digis not produced in simulation
    if (digi.containsSimulation()){
      // look at the early end
      StrawEnd fend = digi.earlyEnd();
      auto const& step =  digi.strawGasStep(fend);
      art::Ptr<SimParticle> const& sp = step->simParticle();
      const unsigned id = sp->id().asInt();
      ++pmap[id];

      const float z = step->position().z();
      const float p = std::sqrt(step->momentum().mag2());
      const float t = step->time();
      if(_debugLevel > 5) printf(" ID = %2i, z = %7.1f, p = %6.2f, t = %6.1f\n",
                                 (int) id, z, p, t);
      if(info_map.find(id) == info_map.end()) {
        info_map[id] = Sim_t();
        info_map[id].id_ = id;
        info_map[id].hit_start_z_ =  1.e6;
        info_map[id].hit_end_z_   = -1.e6;
      }
      auto& info = info_map[id];
      ++info.nhits_;
      if(info.hit_start_z_ > z) {
        info.hit_start_z_ = z;
        info.hit_start_x_ = step->position().x();
        info.hit_start_y_ = step->position().y();
        info.hit_start_p_ = p;
        info.hit_start_pt_ = step->momentum().rho();
        info.hit_start_px_ = step->momentum().x();
        info.hit_start_py_ = step->momentum().y();
        info.hit_start_pz_ = step->momentum().z();
        info.hit_start_t_ = t;
      }
      if(info.hit_end_z_ < z) {
        info.hit_end_z_ = z;
        info.hit_end_x_ = step->position().x();
        info.hit_end_y_ = step->position().y();
        info.hit_end_p_ = p;
        info.hit_end_pt_ = step->momentum().rho();
        info.hit_end_px_ = step->momentum().x();
        info.hit_end_py_ = step->momentum().y();
        info.hit_end_pz_ = step->momentum().z();
        info.hit_end_t_ = t;
      }
    }
  }

  // Store the information for the sim particles of interest
  for(auto sim : *simcol) {
    Sim_t sim_t;
    sim_t.id_ = sim.second.id().asInt();
    if(pmap.find(sim_t.id_) != pmap.end()) {
      const auto& info = info_map[sim_t.id_];
      sim_t = info;
    }
    _simInfo[sim_t.id_] = sim_t;
  }
}

DEFINE_ART_CLASS_TOOL(AgnosticHelixFinderDiag)
