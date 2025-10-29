// Tool to analzye/debug Agnostic Helix Finder processing
// Original author: Matthew Stortini

// Framework
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

// Offline
#include "Offline/CalPatRec/inc/AgnosticHelixFinder_types.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

// ROOT
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TObject.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TGFrame.h"
#include "TRootEmbeddedCanvas.h"
#include "TROOT.h"
#include "TGClient.h"
#include "TApplication.h"
#include "TRint.h"

namespace mu2e {

  using namespace AgnosticHelixFinderTypes;

  class AgnosticHelixFinderDiag : public mu2e::ModuleHistToolBase {

  public:
    enum { kNEventHistsSets = 10,
           kNTimeClusterHistsSets = 10,
           kNHelixSeedHistsSets = 10,
           kNLineInfoHistsSets = 10,
           kNLineSegmentHistsSets = 10,
           kNTripletHistsSets = 10
    };

    struct EventHists {
      TH1* nHelices;
      TH1* nTimeClusters;
      TH1* MC_nSim;
    };

    struct TimeClusterHists {
      TH1* nHelicesPerTC;
      TH1* nComboHitsPerTC;
      TH1* nStrawHitsPerTC;
    };

    struct HelixSeedHists {
      TH1* eDepAvg;
      TH1* p;
      TH1* t0;
      TH1* radius;
      TH1* lambda;
      TH1* d0;
      TH1* hits;
      TH1* MC_PDG;
      TH1* MC_p;
      TH1* MC_dp;
      TH1* MC_purity;
      TH1* MC_hitFrac;
      TH1* MC_hits;
    };

    struct LineInfoHists {
      TH1* nHitsRatioPerHS;
    };

    struct LineSegmentHists {
      TH1* chi2dof;
      TH1* maxHitGap;
    };

    struct TripletHists {
      TH1* radius;
      TH1* rmax;
      TH1* d0;
      TH1* dz12;
      TH1* dz23;
      TH1* dz13;
      TH1* d12;
      TH1* d23;
      TH1* d13;
      TH1* dphi12;
      TH1* dphi23;
      TH1* dphi13;

      TH1* MC_PDG_1;
      TH1* MC_PDG_2;
      TH1* MC_PDG_3;
      TH1* MC_p;
      TH1* MC_pt;
      TH1* MC_dr;
    };

    struct Hists {
      EventHists*       _eventHists      [kNEventHistsSets      ];
      TimeClusterHists* _timeClusterHists[kNTimeClusterHistsSets];
      HelixSeedHists*   _helixSeedHists  [kNHelixSeedHistsSets  ];
      LineInfoHists*    _lineInfoHists   [kNLineInfoHistsSets   ];
      LineSegmentHists* _lineSegmentHists[kNLineSegmentHistsSets];
      TripletHists*     _tripletHists    [kNTripletHistsSets    ];
    };

  private:
    //--------------------------------------------------------------------------------------
    // For sim particle info
    struct Sim_t {
      unsigned id_ = 0; // Sim particle unique ID
      unsigned nhits_ = 0; // Tracker info
      float hit_start_z_ = 0.; // Lowest z tracker hit
      float hit_start_p_ = 0.;
      float hit_start_pt_ = 0.;
      float hit_start_px_ = 0.;
      float hit_start_py_ = 0.;
      float hit_start_pz_ = 0.;
      float hit_start_x_ = 0.;
      float hit_start_y_ = 0.;
      float hit_start_t_ = 0.;
      float hit_end_z_ = 0.; // Highest z tracker hit
      float hit_end_p_ = 0.;
      float hit_end_pt_ = 0.;
      float hit_end_px_ = 0.;
      float hit_end_py_ = 0.;
      float hit_end_pz_ = 0.;
      float hit_end_x_ = 0.;
      float hit_end_y_ = 0.;
      float hit_end_t_ = 0.;
      float start_p_ = 0.; // Info at generation
      float start_pz_ = 0.;
      float start_z_ = 0.;
    };

  public:
    AgnosticHelixFinderDiag(const fhicl::Table<mu2e::AgnosticHelixFinderTypes::Config>& config);
    ~AgnosticHelixFinderDiag();

  private:
    int bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir);
    int bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir);
    int bookHelixSeedHistograms(HelixSeedHists* Hist, art::TFileDirectory* Dir);
    int bookLineInfoHistograms(LineInfoHists* Hist, art::TFileDirectory* Dir);
    int bookLineSegmentHistograms(LineSegmentHists* Hist, art::TFileDirectory* Dir);
    int bookTripletHistograms(TripletHists* Hist, art::TFileDirectory* Dir);

    int fillEventHistograms(EventHists* Hist, diagInfo* Data);
    int fillTimeClusterHistograms(TimeClusterHists* Hist, diagInfo* Data, int loopIndex);
    int fillHelixSeedHistograms(HelixSeedHists* Hist, const hsInfo& info);
    int fillLineInfoHistograms(LineInfoHists* Hist, diagInfo* Data, int loopIndex);
    int fillLineSegmentHistograms(LineSegmentHists* Hist, diagInfo* Data, int loopIndex);
    int fillTripletHistograms(TripletHists* Hist, const tripletInfo& info);

    const SimParticle* helixSimMatch(const HelixSeed& hlx, double& pmc, int& mc_hits);
    const SimParticle* circleSimMatch();
    const SimParticle* tripletSimMatch(const triplet& trip);
    void checkSimTriplets();
    float MCHitPurity(const HelixSeed& hlx);
    float MCHitFraction(const HelixSeed& hlx);
    const SimParticle* hitSim(int index);
    void initSimInfo(const SimParticleCollection* simcol, const StrawDigiMCCollection* digcol);
    int MCColor(int pdg, double pmc = 0.) {
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

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override;
    virtual int fillHistograms(void* Data, int Mode = -1) override;

    //-----------------------------------------------------------------------------
    // functions for display mode
    //-----------------------------------------------------------------------------
    void initDisplay        ();
    void plotCircle         (float xc, float yc, float r,
                             int color, int marker_style,
                             std::string title = "",
                             bool draw_center = true);
    void plotHitXY          (const ComboHit& hit, int color);
    void plotHitXY          (const tripletPoint& point, int color);
    void plotHitPhiZ        (double z, double phi, int color, int iline = 0, int marker = -1);
    void plotHitPhiZ        (int index, int color, int iline = 0, int marker = -1);
    void plotLinePhiZ       (const lineInfo& line, int iline = 0);
    void plotTracker        ();
    void plotHelixXY        (const HelixSeed& hseed, const int index = -1);
    void plotHelixPhiZ      (const HelixSeed& hseed, const int index = -1);
    void plotSimPhiZ        (const SimParticle* sim, const int index = -1);
    void plotCircleXY       ();
    void plotTripletXY      (const tripletInfo& info, const int index = -1, bool mc_triplet = false);
    void plotXY             (int stage);
    void plotPhiZ           (int stage);
    void plotSegment        (int option);
    void plotTotal          (int stage);
    void plotLabel          (float left_margin, float right_margin) {
      TLatex* label = new TLatex();
      label->SetNDC();
      label->SetTextFont(43);
      label->SetTextSize(24);
      label->SetTextAlign(11);
      label->DrawLatex(left_margin, 1.01 - right_margin,
                       Form("Event: %i:%i:%u",
                            _data->event->run(),
                            _data->event->subRun(),
                            _data->event->event()));
      _drawn_objects.push_back(label);
    }

    const SimParticle* simByID(const SimParticleCollection* simcol, const unsigned id) {
      return simcol->getOrNull(cet::map_vector_key(id));
    }
    int MCHits(const SimParticle* sim) {
      if(!sim) return -1;
      const unsigned id = sim->id().asInt();
      if(_simInfo.find(id) == _simInfo.end()) {
        if(_debugLevel > 0) printf("AgnosticHelixFinderDiag::%s: Can't find Sim ID %u in the sim map!\n",
                                   __func__, id);
        return -1;
      }
      return _simInfo[id].nhits_;
    }

    void MCCircle(const int sim_id, float& xC, float& yC, float& r, float& lambda, float& phi0) {
      xC = 0.; yC = 0.; r = 0.; lambda = 0.; phi0 = 0.;

      // Get the sim info
      if(_simInfo.find(sim_id) == _simInfo.end()) return;
      auto& info = _simInfo[sim_id];
      const SimParticle* sim = simByID(_simcol, sim_id);
      if(!sim) return;

      // Evaluate the circle center
      int charge = 1;
      const int pdg = sim->pdgId();
      switch(pdg) {
      case 11: case 13: case -211: charge = -1; break;
      }
      const float MeVmm = 10. / 3. / _data->bz0;
      float xC_1, xC_2, yC_1, yC_2, r_1, r_2, l_1, l_2, phi0_1, phi0_2;
      { // Evaluate using the lowest z hit
        const float x   = info.hit_start_x_ ;
        const float y   = info.hit_start_y_ ;
        const float z   = info.hit_start_z_ ;
        const float t   = info.hit_start_t_ ;
        const float px  = info.hit_start_px_;
        const float py  = info.hit_start_py_;
        const float pz  = info.hit_start_pz_;
        const float pt  = info.hit_start_pt_;
        const float phi = polyAtan2(px, py) + charge*M_PI/2.;
        r_1  = pt * MeVmm; // convert to mm
        xC_1 = x + charge*r_1*py/pt; //std::sin(phi);
        yC_1 = y - charge*r_1*px/pt; //std::cos(phi);
        l_1 = pz * MeVmm;
        // phi = phi0 + z / lambda
        phi0_1 = std::fmod(phi - z / l_1, M_PI);
        if(_debugLevel > 2) printf(" MC info: Start: x = (%6.1f, %6.1f, %7.1f, %7.1f), p = (%6.1f, %6.1f, %6.1f), pt = %5.1f --> r = %5.1f, center = (%6.1f, %6.1f), lambda = %6.1f, phi0 = %4.1f\n",
                                   x, y, z, t, px, py, pz, pt, r_1, xC_1, yC_1, l_1, phi0_1);
      }
      { // Evaluate using the highest z hit
        const float x   = info.hit_end_x_ ;
        const float y   = info.hit_end_y_ ;
        const float z   = info.hit_end_z_ ;
        const float t   = info.hit_end_t_ ;
        const float px  = info.hit_end_px_;
        const float py  = info.hit_end_py_;
        const float pz  = info.hit_end_pz_;
        const float pt  = info.hit_end_pt_;
        const float phi = polyAtan2(px, py) + charge*M_PI/2.;
        r_2  = pt * MeVmm; // convert to mm
        xC_2 = x + charge*r_2*py/pt;//std::sin(phi);
        yC_2 = y - charge*r_2*px/pt;//std::cos(phi);
        l_2 = pz * MeVmm;
        phi0_2 = std::fmod(phi - z / l_2, M_PI);
        if(_debugLevel > 2) printf(" MC info: End: x = (%6.1f, %6.1f, %7.1f, %7.1f), p = (%6.1f, %6.1f, %6.1f), pt = %5.1f --> r = %5.1f, center = (%6.1f, %6.1f), lambda = %6.1f, phi0 = %4.1f\n",
                                   x, y, z, t, px, py, pz, pt, r_2, xC_2, yC_2, l_2, phi0_2);
      }
      // Just take the average of the two point estimates
      xC     = 0.5*( xC_1 +  xC_2);
      yC     = 0.5*( yC_1 +  yC_2);
      r      = 0.5*(  r_1 +   r_2);
      lambda = 0.5*(  l_1 +   l_2);
      phi0   = phi0_1; // FIXME: Add logic for angle bisector
    }
    void MCCircle(const int sim_id, float& xC, float& yC, float& r) {
      float lambda, phi0;
      MCCircle(sim_id, xC, yC, r, lambda, phi0);
    }

    const SimParticle* helixSimMatch(const HelixSeed& hlx, double& pmc) {
      int mc_hits;
      return helixSimMatch(hlx, pmc, mc_hits);
    }
    const SimParticle* helixSimMatch(const HelixSeed& hlx, int& mc_hits) {
      double pmc;
      return helixSimMatch(hlx, pmc, mc_hits);
    }
    const SimParticle* helixSimMatch(const HelixSeed& hlx) {
      int mc_hits; double pmc;
      return helixSimMatch(hlx, pmc, mc_hits);
    }

    int fillHelixSeedHistograms(HelixSeedHists* Hist, diagInfo* Data, int loopIndex) {
      if(loopIndex >= int(Data->helixSeedData.size())) return -1;
      const auto& info = Data->helixSeedData.at(loopIndex);
      return fillHelixSeedHistograms(Hist, info);
    }

  protected:
    std::string _simTag;
    std::string _digTag;
    bool        _display;
    bool        _display3D;
    bool        _showProtons;
    double      _pMin;
    double      _tMin;
    int         _debugLevel;

    Hists _hist;
    diagInfo* _data;
    const SimParticleCollection* _simcol = nullptr;
    const StrawDigiMCCollection* _digcol = nullptr;

    std::map<unsigned, Sim_t> _simInfo;
    std::map<unsigned, tripletInfo> _simTriplets; // map sims to their closest triplet

    // For running a display
    TCanvas* c_tot_ = nullptr;
    TCanvas* c_3d_  = nullptr;
    TPad* c_hlx_ = nullptr;
    TPad* c_seg_ = nullptr;
    TPad* c_trp_ = nullptr;
    std::vector<TObject*> _drawn_objects;
    TGMainFrame* _gMainFrame = nullptr;
    TApplication* _application = nullptr;

  };

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
      if(_display) initDisplay();
      return 0;
    }

    //-----------------------------------------------------------------------------
    // Triplet stage histograms
    //-----------------------------------------------------------------------------

    if(Mode == DIAG::kTriplet) {
      if(_display && _data->diagLevel > 5) {
        std::cout << "Triplet stage: Loop condition " << ConditionName(_data->loopCondition)
                  << std::endl;
        plotXY(0);
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
        plotXY(2);
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
        plotXY(1); // triplet
        plotSegment(0); // line segments
      }
      fillTripletHistograms(_hist._tripletHists[1], _data->tripInfo);
      return 0;
    }

    //-----------------------------------------------------------------------------
    // Line stage histograms
    //-----------------------------------------------------------------------------

    if(Mode == DIAG::kLine) {
      if(_display && _data->diagLevel > 3) {
        std::cout << "Line stage: Loop condition " << ConditionName(_data->loopCondition)
                  << std::endl;
        plotXY(1); // triplet
        plotSegment(0); // line segments
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
        if(_display && _data->diagLevel > 2) {
          std::cout << "Helix stage: Loop condition " << ConditionName(_data->loopCondition)
                    << std::endl;
          plotXY(1); // plot the triplet
          plotSegment(2); // plot the helix phi z
          plotXY(3); // plot the helix x y
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
        plotSegment(3);
        plotXY(5);
        plotXY(6); // plot the helices in the time cluster
      }
      return 0;
    }

    //-----------------------------------------------------------------------------
    // Final histograms, only filled at event end
    //-----------------------------------------------------------------------------

    if(Mode == DIAG::kEnd) {
      if(_display) {
        plotSegment(3);
        plotXY(5);
        plotXY(4);
        if(_display3D) plotTotal(0);
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
  const SimParticle* AgnosticHelixFinderDiag::helixSimMatch(const HelixSeed& hlx, double& pmc, int& mc_hits) {
    pmc = -1.; mc_hits = -1;

    // Ensure valid collections
    if(hlx.hits().empty()) return nullptr;
    if(!_simcol) return nullptr;
    if(!_digcol) return nullptr;

    // Loop through the hits on the helix looking for the sim with the most hits
    const auto& hits = hlx.hits();
    const size_t nhits = hits.size();
    if(nhits == 0) return nullptr;
    int max_hits(-1), max_id(-1); // current maxima
    std::map<int, int> id_to_hits;
    for(size_t ihit = 0; ihit < nhits; ++ihit) {
      std::vector<StrawDigiIndex> shids;
      hits.fillStrawDigiIndices(ihit,shids);
      for (size_t j = 0; j < shids.size(); ++j) {
        if(shids[j] >= _digcol->size()) continue;
        const auto& sdmc   = _digcol->at(shids[j]);
        const auto& spmcp  = sdmc.earlyStrawGasStep();
        const auto& simptr = spmcp->simParticle();
        const int sim_id   = simptr->id().asInt();
        ++id_to_hits[sim_id];
        if(id_to_hits[sim_id] > max_hits) {
          max_hits = id_to_hits[sim_id];
          max_id = sim_id;
          mc_hits = max_hits;
          pmc = std::sqrt(spmcp->momentum().mag2());
        }
      }
    }

    // Return the sim with the most hits
    const SimParticle* sim = (max_id >= 0) ? simByID(_simcol, max_id) : nullptr;
    return sim;
  }

  //--------------------------------------------------------------------------------------
  const SimParticle* AgnosticHelixFinderDiag::circleSimMatch() {
    // Ensure valid collections
    if(!_data || !_data->tcHits || !_data->chColl) return nullptr;
    if(!_simcol) return nullptr;
    if(!_digcol) return nullptr;

    // Loop through the hits on the circle looking for the sim with the most hits
    int max_hits(-1), max_id(-1); // current maxima
    std::map<int, int> id_to_hits;
    for(size_t ihit = 0; ihit < _data->tcHits->size(); ++ihit) {
      const auto& tcHit = _data->tcHits->at(ihit);
      if(tcHit.inHelix || !tcHit.used || tcHit.hitIndice < 0) continue;
      std::vector<StrawDigiIndex> shids;
      _data->chColl->fillStrawDigiIndices(tcHit.hitIndice,shids);
      for (size_t j = 0; j < shids.size(); ++j) {
        if(shids[j] >= _digcol->size()) continue;
        const auto& sdmc   = _digcol->at(shids[j]);
        const auto& spmcp  = sdmc.earlyStrawGasStep();
        const auto& simptr = spmcp->simParticle();
        const int sim_id   = simptr->id().asInt();
        ++id_to_hits[sim_id];
        if(id_to_hits[sim_id] > max_hits) {
          max_hits = id_to_hits[sim_id];
          max_id = sim_id;
        }
      }
    }

    // Return the sim with the most hits
    const SimParticle* sim = (max_id >= 0) ? simByID(_simcol, max_id) : nullptr;
    return sim;
  }

  //--------------------------------------------------------------------------------------
  const SimParticle* AgnosticHelixFinderDiag::tripletSimMatch(const triplet& trip) {

    // Sim particle for each hit
    const SimParticle* sim_1 = hitSim(trip.i.hitIndice);
    const SimParticle* sim_2 = hitSim(trip.j.hitIndice);
    const SimParticle* sim_3 = hitSim(trip.k.hitIndice);

    // Sim IDs to see if one sim made multiple hits
    const int sim_id_1 = (sim_1) ? sim_1->id().asInt() : -1;
    const int sim_id_2 = (sim_2) ? sim_2->id().asInt() : -1;
    const int sim_id_3 = (sim_3) ? sim_3->id().asInt() : -1;

    const SimParticle* best_sim = nullptr;
    if     (sim_1 && sim_2 && sim_id_1 == sim_id_2) best_sim = sim_1; // Take the sim with more hits
    else if(sim_1 && sim_3 && sim_id_1 == sim_id_3) best_sim = sim_1;
    else if(sim_2 && sim_3 && sim_id_2 == sim_id_3) best_sim = sim_2;
    else if(sim_1)                                  best_sim = sim_1; // None match, take an existing one
    else if(sim_2)                                  best_sim = sim_2;
    else if(sim_3)                                  best_sim = sim_3;

    return best_sim;
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
  const SimParticle* AgnosticHelixFinderDiag::hitSim(int index) {
    if(!_simcol || !_digcol || !_data || !_data->chColl) return nullptr; // missing necessary collections

    // Get the MC info for the hit
    if(index < 0) return nullptr; // calo cluster/stopping target use non-physical indices

    // Get the sim with the most hits in this combo hit
    std::vector<StrawDigiIndex> shids;
    _data->chColl->fillStrawDigiIndices(index, shids);
    int max_hits(-1), max_id(-1);
    std::map<int, int> id_to_hits;
    for (size_t i = 0; i < shids.size(); ++i) {
      if(shids[i] >= _digcol->size()) continue;
      const auto& sdmc   = _digcol->at(shids[i]);
      const auto& spmcp  = sdmc.earlyStrawGasStep();
      const auto& simptr = spmcp->simParticle();
      const int sim_id   = simptr->id().asInt();
      ++id_to_hits[sim_id];
      if(id_to_hits[sim_id] > max_hits) {
        max_hits = id_to_hits[sim_id];
        max_id = sim_id;
      }
    }
    return simByID(_simcol, max_id);
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

  //-----------------------------------------------------------------------------
  void AgnosticHelixFinderDiag::initDisplay() {

    // Need the TApplication for the display, where TRint allows command-line usage
    int    tmp_argc(0);
    char** tmp_argv(0);
    if(!gApplication && !_application) {
      _application = new TRint("AgnosticHelixFinderDiag_tool", &tmp_argc, tmp_argv);
      gApplication = _application;
    } else if(gApplication && !_application) _application = gApplication; // use the existing one

    // Ensure the display canvases are available
    if(_display3D && !c_3d_) {
      c_3d_ = new TCanvas("c_3d", "3D canvas"  , 900, 900);
    }
    if(!c_tot_) {
      c_tot_ = new TCanvas("c_tot", "Total canvas"  , 1260, 900);
    }
    if(!c_hlx_) {
      c_hlx_ = new TPad("c_hlx", "Helix canvas"  , 0., 0.3, 0.5, 1.);
      c_hlx_->SetLeftMargin(0.12);
      c_hlx_->SetBottomMargin(0.10);
      c_hlx_->SetTicks(1);
      c_tot_->cd();
      c_hlx_->Draw();
    }
    if(!c_seg_) {
      c_seg_ = new TPad("c_seg", "Segment canvas"  , 0., 0.0, 1.0, 0.3);
      c_seg_->SetLeftMargin(0.06);
      c_seg_->SetBottomMargin(0.22);
      c_seg_->SetRightMargin(0.01);
      c_seg_->SetTopMargin(0.08);
      c_seg_->SetTicks(1);
      c_tot_->cd();
      c_seg_->Draw();
    }
    if(!c_trp_) {
      c_trp_ = new TPad("c_trp", "Triplet canvas"  , 0.5, 0.3, 1.0, 1.0);
      c_trp_->SetLeftMargin(0.12);
      c_trp_->SetBottomMargin(0.12);
      c_trp_->SetTicks(1);
      c_tot_->cd();
      c_trp_->Draw();
    }

    // clean up the last event
    for(auto o : _drawn_objects) delete o;
    _drawn_objects.clear();
  }

  //-----------------------------------------------------------------------------
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
  void AgnosticHelixFinderDiag::plotLinePhiZ(const lineInfo& line, int iline) {

    constexpr int colors[] = {kBlue, kGreen, kRed, kOrange, kPink, kYellow-2};
    const static int ncolors = sizeof(colors) / sizeof(*colors);
    const int color = colors[iline % ncolors];

    // Plot the points
    for(auto index : line.tcHitsIndices) plotHitPhiZ(index, color, iline);

    // Draw the fit line
    const float zmin  = line.zMin;
    const float zmax  = line.zMax;
    // clone the fitter to avoid changing its state
    ::LsqSums2 fitter = line.fitter;
    const float slope = fitter.dydx();
    const float phi0  = fitter.y0();;
    const float phi1  = slope * zmin + phi0;
    const float phi2  = slope * zmax + phi0;
    TLine* tline = new TLine(zmin, phi1, zmax, phi2);
    tline->SetLineColor(color);
    tline->SetLineWidth(1);
    tline->Draw("same");

    _drawn_objects.push_back(tline);
  }

  //-----------------------------------------------------------------------------
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
        const int mc_color = (pdg == 2212) ? kOrange : kBlue;
        plotCircle(xC_mc, yC_mc, r_mc, mc_color, 2, title); // MC helix estimate
      }
    }
  }

  //-----------------------------------------------------------------------------
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
  void AgnosticHelixFinderDiag::plotSimPhiZ(const SimParticle* sim, const int index) {
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
        plotHitPhiZ(z, phi0 + z/lambda, color, iline, 2);
      }
    }
  }

  //-----------------------------------------------------------------------------
  void AgnosticHelixFinderDiag::plotCircleXY() {
    if(!_data || !_data->tcHits || !_data->chColl) return;

    // check for MC info, plot below the reco data
    const auto sim = circleSimMatch();
    const int sim_id = (sim) ? sim->id().asInt() : -1;
    if(sim && _simInfo.find(sim_id) != _simInfo.end()) {
      float xC_mc, yC_mc, r_mc;
      MCCircle(sim_id, xC_mc, yC_mc, r_mc);
      const int pdg = std::abs(sim->pdgId());
      const int mc_color = (pdg == 2212) ? kOrange : kBlue;
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
  void AgnosticHelixFinderDiag::plotTripletXY(const tripletInfo& info, const int index, bool mc_triplet) {

    // Plot the circle from the triplet
    const float r  = info.radius;
    const float xC = info.xC;
    const float yC = info.yC;
    if(_debugLevel > 2) printf("    Triplet   %i: x = %6.1f, y = %6.1f, r = %6.1f\n", index, xC, yC, r);
    std::string title = (index >= 0) ? Form("TripletReco_%i", index) : "TripletReco";
    plotCircle(xC, yC, r, kRed, 2, title); // Reco triplet

    // Plot the triplet points used
    const int color = (mc_triplet) ? kOrange : kGreen;
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
  // function to plot XY view at various stages of logic
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinderDiag::plotXY(int stage) {
    if(_debugLevel > -1) printf("  %12s: Plotting XY for stage %i\n", __func__, stage);
    if(!_data) return;
    if(!c_tot_) return;

    TPad* c(nullptr);
    switch(stage) {
    case 0: c = c_trp_; break;
    case 1: c = c_trp_; break; // triplet, but don't stop to animate
    case 5: c = c_trp_; break;
    case 2: c = c_hlx_; break; // seed circle
    case 3: c = c_hlx_; break; // final helix
    case 4: c = c_hlx_; break; // total event
    case 6: c = c_hlx_; break; // time cluster stage
    }
    if(!c) {
      printf("[AgnosticHelixFinderDiag::%s] Undefined fit stage %i\n", __func__, stage);
      return;
    }
    c_tot_->cd();
    c->Draw();
    c->cd();

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

    // Stage 0/1: Per triplet
    if(stage == 0 || stage == 1) plotTripletXY(_data->tripInfo);

    // Stage 2: Circle fit
    if(stage == 2) plotCircleXY();

    // Stage 3: Per helix
    if(stage == 3 && _data->hseed) plotHelixXY(*(_data->hseed));

    // Stage 4: End of the event processing: helix canvas
    if(stage == 4 || stage == 6) {

      // All hits in the event
      if(_data && _data->chColl && (stage != 6 || _data->tc)) {
        const size_t nhits = (stage == 6) ? _data->tc->hits().size() : _data->chColl->size();
        for(size_t ihit = 0; ihit < nhits; ++ihit) {
          const size_t index = (stage == 6) ? _data->tc->hits().at(ihit) : ihit;
          const auto& hit = _data->chColl->at(index);
          if(hit.time() < _tMin) continue;
          const auto sim = hitSim(index);
          double pmc(0.);
          int pdg = 0;
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

      // All helices
      int index = 0;
      for(auto& info : _data->helixSeedData) {
        plotHelixXY(info.seed, index);
        ++index;
      }
    }

    // Stage 4/5: End of the event processing: MC truth
    if(stage == 4 || stage == 5) {
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
            if(stage == 4) { // helix canvas
              float xC_mc, yC_mc, r_mc;
              MCCircle(info.id_, xC_mc, yC_mc, r_mc);
              if(_debugLevel > 2) printf(" Relevant MC: pdg = %4i, pmc = %6.1f, x(C) = %7.1f, y(C) = %7.1f, r = %6.1f\n",
                                         pdg, pmc, xC_mc, yC_mc, r_mc);
              const int mc_color = MCColor(pdg, pmc);
              plotCircle(xC_mc, yC_mc, r_mc, mc_color, 2, "MC_Sim"); // MC helix estimate
            } else if(stage == 5) { // triplet canvas
              if(_simTriplets.find(info.id_) != _simTriplets.end()) {
                plotTripletXY(_simTriplets[info.id_], index, true);
                ++index;
              }
            }
          }
        }
      }
    }

    plotLabel(c->GetLeftMargin(), c->GetTopMargin());
    c->Modified(); c->Update();
    if(stage == 0 || stage == 3 || (!_display3D && stage == 4) || stage == 6) {
      _application->Run(true);
      std::cout << std::endl;
    }

    if(_debugLevel > 3) c->Print(Form("xy_display_%i.png", stage));
    if(_debugLevel > 1) printf("  %s: Finished plotting XY for stage %i\n", __func__, stage);
  }

  //-----------------------------------------------------------------------------
  // function to plot a phi-z segment
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinderDiag::plotSegment(int option) {
    if(_debugLevel > -1) printf("  %12s: Plotting segment with option %i\n", __func__, option);
    if(!_data || !_data->tcHits) return;
    if(!c_tot_) return;

    // Retrieve the relevant canvas and switch to this pad
    TPad* c(nullptr);
    switch(option) {
    case 0: c = c_seg_; break;
    case 1: c = c_seg_; break;
    case 2: c = c_seg_; break; // plot helix PhiZ
    case 3: c = c_seg_; break; // all helices PhiZ
    case 4: c = c_seg_; break; // all time cluster helices PhiZ
    }
    if(!c) {
      printf("AgnosticHelixFinderDiag::%s: Undefined option %i\n", __func__, option);
      return;
    }
    c_tot_->cd();
    c->Draw();
    c->cd();

    // Use a graph to control the axes
    TGraph* axes = new TGraph();
    axes->AddPoint(0.,0.); // necessary for drawing
    axes->SetMarkerSize(0.);
    axes->GetXaxis()->SetLimits(-1600., 1600.);
    axes->GetYaxis()->SetRangeUser(-7.0, 7.0);
    axes->GetXaxis()->SetTitle("z (mm)");
    axes->GetYaxis()->SetTitle("#phi");
    axes->Draw("AP");
    _drawn_objects.push_back(axes);

    axes->GetXaxis()->SetLabelSize(0.08);
    axes->GetYaxis()->SetLabelSize(0.08);
    axes->GetXaxis()->SetTitleSize(0.10);
    axes->GetYaxis()->SetTitleSize(0.10);
    axes->GetYaxis()->SetTitleOffset(0.3);

    // Plot each line information
    if(option == 0 || option == 1) {
      if(!_data->seedPhiLines) {
        printf("[AgnosticHelixFinderDiag::%s] Exiting due to missing data!\n", __func__);
        return;
      }
      for(size_t iline = 0; iline < _data->seedPhiLines->size(); ++iline) {
        plotLinePhiZ(_data->seedPhiLines->at(iline), iline);
      }
    } else if(option == 2 || option == 3) { // helix info
      if(option == 2 && _data->hseed) plotHelixPhiZ(*(_data->hseed));
      if(option == 3) {
        // Relevant sim particles
        int index = 0;
        for(auto info_pair : _simInfo) {
          auto& info = info_pair.second;
          const double pmc = std::max(info.hit_start_p_, info.hit_end_p_);
          const double tmc = std::max(info.hit_start_t_, info.hit_end_t_);
          if(info.nhits_ > 10 && pmc > _pMin && tmc > _tMin) {
            const auto sim = simByID(_simcol, info.id_);
            plotSimPhiZ(sim, index); // MC helix estimate
            ++index;
          }
        }

        // Helices
        index = 0;
        for(auto& info : _data->helixSeedData) {
          plotHelixPhiZ(info.seed, index);
          ++index;
        }
      }
    }

    plotLabel(c->GetLeftMargin(), c->GetTopMargin());
    c->Modified(); c->Update();
    if(option == 0 || option == 1) _application->Run(true);
  }

  //-----------------------------------------------------------------------------
  // function to plot the 3D event
  //-----------------------------------------------------------------------------
  void AgnosticHelixFinderDiag::plotTotal(int option) {
    if(_debugLevel > -1) printf("  %12s: Plotting with option %i\n", __func__, option);
    if(!_data || !_data->tcHits || !_data->chColl) return;

    // Retrieve the relevant canvas and switch to this pad
    TCanvas* c(nullptr);
    switch(option) {
    case 0: c = c_3d_; break;
    }
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

    // FIXME: Currently 3D drawing throws warnings due to ROOT issue and art upgrades to an error
    if(false) {
      frame3d->Draw();
      g_tracker->Draw("P SAME");
      g_hits->Draw("P SAME");
    } else { // only draw the hits for now
      g_hits->Draw("P SAME");
    }

    plotLabel(c->GetLeftMargin(), c->GetTopMargin());
    c->Modified(); c->Update();
    _application->Run(true);
  }

  DEFINE_ART_CLASS_TOOL(AgnosticHelixFinderDiag)

} // namespace mu2e
