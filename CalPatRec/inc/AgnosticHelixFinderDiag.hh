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

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override;
    virtual int fillHistograms(void* Data, int Mode = -1) override;

    // MC utilities
    void checkSimTriplets();
    float MCHitPurity(const HelixSeed& hlx);
    float MCHitFraction(const HelixSeed& hlx);
    const SimParticle* helixSimMatch(const HelixSeed& hlx, double& pmc, int& mc_hits);
    const SimParticle* circleSimMatch();
    const SimParticle* tripletSimMatch(const triplet& trip);
    const SimParticle* hitSim(int index);
    int MCHits(const SimParticle* sim);
    void MCCircle(const int sim_id, float& xC, float& yC, float& r, float& lambda, float& phi0);
    void initSimInfo(const SimParticleCollection* simcol, const StrawDigiMCCollection* digcol);
    const SimParticle* simByID(const SimParticleCollection* simcol, const unsigned id) {
      return simcol->getOrNull(cet::map_vector_key(id));
    }

    //-----------------------------------------------------------------------------
    // functions for display mode
    //-----------------------------------------------------------------------------
    int MCColor             (int pdg, double pmc = 0.);
    void initDisplay        ();
    void plotCircle         (float xc, float yc, float r,
                             int color, int marker_style,
                             std::string title = "",
                             bool draw_center = true);
    void plotHitXY          (const ComboHit& hit, int color);
    void plotHitXY          (const tripletPoint& point, int color);
    void plotHitPhiZ        (double z, double phi, int color, int iline = 0, int marker = -1);
    void plotHitPhiZ        (int index, int color, int iline = 0, int marker = -1);
    void plotLinePhiZ       (const lineInfo& line, bool resolve = false, int iline = 0);
    void plotTracker        ();
    void plotHelixXY        (const HelixSeed& hseed, const int index = -1);
    void plotHelixPhiZ      (const HelixSeed& hseed, const int index = -1);
    void plotSimPhiZ        (const SimParticle* sim, bool resolve = false, const int index = -1);
    void plotCircleXY       ();
    void plotCirclePhiZ     ();
    void plotTripletXY      (const tripletInfo& info, const int index = -1, bool mc_triplet = false);
    void plotHitsXY         (const bool use_tc);
    void plotMC             (int stage, bool phiz = false);
    void plotBeginStage     (const bool use_tc = false);
    void plotTripletStage   (const bool mc_triplets = false);
    void plotCircleStage    ();
    void plotSegmentStage   (bool resolve = false, bool seed_cirlce = false);
    void plotHelixStageXY   (int stage);
    void plotHelixStagePhiZ (int stage);
    void plotXYAxes         ();
    void plotPhiZAxes       (double phi_min, double phi_max);
    void plotTotal          (int stage);
    void plotLabel          (float left_margin, float top_margin,
                             float right_margin = 0.1, std::string title = "");
    void runDisplay() {
      if(_application) {
        std::cout << "Enter \".q\" to continue or \".qqqq\" to quit\n";
        _application->Run(true); std::cout << std::endl;
      }
    }
    bool beginPlot(TPad* c) {
      if(!_data || !c_tot_ || !c) return false;
      c_tot_->cd(); c->Draw(); c->cd();
      return true;
    }
    void endPlot(TPad* c, std::string title = "") {
      if(!c) return;
      plotLabel(c->GetLeftMargin(), c->GetTopMargin(), c->GetRightMargin(), title);
      c->Modified(); c->Update();
    }


    // Overloaded function calls
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

  }; // end class

  //--------------------------------------------------------------------------------------
  int AgnosticHelixFinderDiag::MCHits(const SimParticle* sim) {
    if(!sim) return -1;
    const unsigned id = sim->id().asInt();
    if(_simInfo.find(id) == _simInfo.end()) {
      if(_debugLevel > 0) printf("AgnosticHelixFinderDiag::%s: Can't find Sim ID %u in the sim map!\n",
                                 __func__, id);
      return -1;
    }
    return _simInfo[id].nhits_;
  }

  //--------------------------------------------------------------------------------------
  void AgnosticHelixFinderDiag::MCCircle(const int sim_id, float& xC, float& yC,
                                         float& r, float& lambda, float& phi0) {
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
  const SimParticle* AgnosticHelixFinderDiag::hitSim(int index) {
    if(!_simcol || !_digcol || !_data || !_data->chColl) return nullptr; // missing necessary collections

    // Get the MC info for the hit
    if(index < 0) return nullptr; // calo cluster/stopping target use non-physical indices

    // Get the sim with the most hits in this combo hit
    std::vector<StrawDigiIndex> shids;
    _data->chColl->fillStrawDigiIndices(index, shids);
    int max_hits(-1), max_id(-1);
    float max_pmc(0.f); // if equal hits, take the higher momentum sim
    std::map<int, int> id_to_hits;
    for (size_t i = 0; i < shids.size(); ++i) {
      if(shids[i] >= _digcol->size()) continue;
      const auto& sdmc   = _digcol->at(shids[i]);
      const auto& spmcp  = sdmc.earlyStrawGasStep();
      const auto& simptr = spmcp->simParticle();
      const int sim_id   = simptr->id().asInt();
      const float pmc = std::sqrt(spmcp->momentum().mag2());
      ++id_to_hits[sim_id];
      if(id_to_hits[sim_id] > max_hits || (id_to_hits[sim_id] == max_hits && pmc > max_pmc)) {
        max_hits = id_to_hits[sim_id];
        max_id = sim_id;
        max_pmc = pmc;
      }
    }
    return simByID(_simcol, max_id);
  }

} // end namespace mu2e
