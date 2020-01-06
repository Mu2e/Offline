//Author: S Middleton April 2019
//Purpose : to help select muons for cosmic study based on momentum and other info.
//This also provides  the MC truth information for the cosmic muon diagnostics

#include "DataProducts/inc/PDGCode.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TNtuple.h"
#include "Math/VectorUtil.h"
using namespace std;

namespace mu2e {

  class CosmicMuonInfo : public art::EDFilter {
  public:

    explicit CosmicMuonInfo(fhicl::ParameterSet const& pset);
    virtual bool filter  (art::Event& event) override;

  private:
    //Config _conf;
    art::InputTag _strawHitsTag, _strawDigisTag, _strawDigiMCsTag, _caloDigisTag;
    art::InputTag _panelHitsTag;
    int           _diagLevel;

    struct Cuts {
      double   pmin;
      double   pmax;
      unsigned minStrawDigis;
      unsigned minPlanes;
      int      minBackground;
      int      maxBackground;
      Cuts ( fhicl::ParameterSet const& pset):
        pmin(pset.get<double>("pmin")),
	pmax(pset.get<double>("pmax")),
        minStrawDigis(pset.get<unsigned>("minStrawDigis")),
        minPlanes(pset.get<unsigned>("minPlanes")),
        minBackground(pset.get<int>("minBackground")),
        maxBackground(pset.get<int>("maxBackground")){
      }
    };
    Cuts _cuts;
    TH1D* _phiMC          = nullptr;
    TH1D* _thetaMC       = nullptr;
    TH1D* _phiMCcuts        = nullptr;
    TH1D* _thetaMCcuts       = nullptr;
    TH1D* _hNStrawHits    = nullptr;
    TH1D* _hNPanelHits    = nullptr;
    TH1D* _hUniquePanel   = nullptr;
    TH1D* _hPlane         = nullptr;
    TH1D* _hPanel         = nullptr;
    TH1D* _hViewId        = nullptr;
    TH1D* _hStraw         = nullptr;
    TH1D* _hnBadDigis     = nullptr;
    TH1D* _hnComboStraw   = nullptr;
    TH1D* _hnComboPanel   = nullptr;
    TH1D* _hnMCTracks     = nullptr;
    TH1D* _hnMCMuons      = nullptr;
    TH1D* _hnDigisPerMuon = nullptr;
    TH1D* _hMomentumAll   = nullptr;
    TH1D* _hnBackground   = nullptr;
    TH1D* _hMomentum      = nullptr;
    TH1D* _hMomentumStart = nullptr;
    TH1D* _hMomentumDelta = nullptr;
    TH1D* _hnPlanes       = nullptr;
    TH1D* _hWireResStraw  = nullptr;
    TH1D* _hWireResPanel  = nullptr;
    TH1D* _hTransResStraw = nullptr;
    TH1D* _hTransResPanel = nullptr;
    TH1D* _hUniquePanelGoodMu   = nullptr;
    TH1D* _hPlaneGoodMu         = nullptr;
    TH1D* _hPanelGoodMu         = nullptr;
    TH1D* _hViewIdGoodMu        = nullptr;
    TH1D* _hStrawGoodMu         = nullptr;
    TH1D* _hFracPanelHits       = nullptr;
    TH1D* _hnMissingPanelHits   = nullptr;

    TH1D* _hnCaloDigiAll  = nullptr;
    TH1D* _hnCaloDigi     = nullptr;

    TNtuple* _ntWireRes   = nullptr;

  };
}

namespace {
  class DigisBySim {

  public:

    // Some MC truth information about tracks that make digis.
    struct TrkMCInfo{

      // The int inside the vector is the common index into these collections:
      // StrawDigi, StrawHit, StrawDigiMC and the ComboHitCollection made by makeSH.
      std::vector<int> digi_indices;

      // Momentum of the track; max of MC truth momentum over all digis.
      double p;

      TrkMCInfo():digi_indices(),p(0){}

      void addIndex( size_t i, double pStep ){
        p = std::max(p,pStep);
        digi_indices.push_back(i);
      }

      int at( size_t i ) const{
        return digi_indices.at(i);
      }

      size_t size() const{
        return digi_indices.size();
      }

    };
    typedef std::map<art::Ptr<mu2e::SimParticle>,TrkMCInfo > impl_type;

    DigisBySim( mu2e::StrawDigiMCCollection const& digimcs ):
      _bySim(){

      int n{-1};
      for ( auto const& digimc : digimcs ){
        auto const& sim_cal = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
        auto const& sim_hv  = digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle();
        if ( sim_cal == sim_hv ){
          double p = sqrt(digimc.strawGasStep(mu2e::StrawEnd::cal)->momentum().mag2());
          _bySim[digimc.strawGasStep(mu2e::StrawEnd::cal)->simParticle()].addIndex(++n,p);
          ++_sum;
        }else{
          ++_nbad;
        }
      }
    }

    impl_type const& bySim() const { return _bySim; }
    int sum()  const { return _sum; }
    int nBad() const { return _nbad; }
    int nTracks() const  { return _bySim.size(); }

    int maxDigis() const {
      int maxD{0};
      for ( auto const& i : _bySim ){
        maxD = std::max( maxD, int(i.second.size()) );
      }
      return maxD;
    }

    int nGoodTracks ( int minDigis ){
      int ngood{0};
      for ( auto const& i : _bySim ){
        if ( int(i.second.digi_indices.size()) >= minDigis ){
          ++ngood;
        }
      }
      return ngood;
    }

    int nMuons () const{
      int nMu{0};
      for ( auto const& i : _bySim ){
        mu2e::SimParticle const& sim = *i.first;
        if ( std::abs(sim.pdgId()) == mu2e::PDGCode::mu_minus ){
          ++nMu;
        }
      }
      return nMu;
    }

  private:
    impl_type  _bySim;
    int        _sum  = 0; // Sum of all digis over all SimParticles
    int        _nbad = 0; // Number of digis rejected since opposite ends
    // of the wire were set by different particles.
  };


  // Define an id for each of the 12 views.
  // Plane rotation pattern is AB BA AB BA ...
  int viewId( mu2e::StrawId id){
    int modStation = id.getStation()%2;
    int modPlane   = id.getPlane()%2;
    int swapPlane  = (modPlane == 0) ? 1 : 0;
    if( modStation == 0 ){
      return modPlane*mu2e::StrawId::_npanels + id.getPanel();
    }
    return swapPlane*mu2e::StrawId::_npanels + id.getPanel();
  }

}

mu2e::CosmicMuonInfo::CosmicMuonInfo(fhicl::ParameterSet const& pset):
  art::EDFilter{pset},
  _strawHitsTag(pset.get<std::string>("strawHitsTag")),
  _strawDigisTag(pset.get<std::string>("strawDigisTag")),
  _strawDigiMCsTag(pset.get<std::string>("strawDigiMCsTag")),
  _caloDigisTag(pset.get<std::string>("caloDigisTag")),
  _panelHitsTag(pset.get<std::string>("panelHitsTag")),
  _diagLevel(pset.get<int>("diagLevel",0)),
  _cuts(pset.get<fhicl::ParameterSet>("filterCuts")){

  art::ServiceHandle<art::TFileService> tfs;
  
  _phiMC	  = tfs->make<TH1D>( "#phi_{MC}",   "Angle #phi_{MC} of Muon MC tracks All",  100, -3.141529,      3.141529 );
  _thetaMC        = tfs->make<TH1D>( "#theta_{MC}",   "Angle #theta_{MC} of Muon MC tracks All",  20, 0,      3.141529 );
  _phiMCcuts	  = tfs->make<TH1D>( "#phi_after_cuts_{MC}",   "Angle #phi_{MC} of Muon MC tracks after MC cuts",  100, -3.141529,      3.141529 );
  _thetaMCcuts        = tfs->make<TH1D>( "#theta_after_cuts{MC}",   "Angle #theta_{MC} of Muon MC tracks aft MC cuts",  20, 0,      3.141529 );
  _hNStrawHits    = tfs->make<TH1D>( "hNStrawHits",   "Number of Straw Hits",         100,  0.,       100. );
  _hNPanelHits    = tfs->make<TH1D>( "hNPanelHits",   "Number of Panel Hits",         100,  0.,       100. );
  _hUniquePanel   = tfs->make<TH1D>( "hUniquePanel",  "Unique Panel ID",              216,  0.,       216. );
  _hPlane         = tfs->make<TH1D>( "hPlane",        "Plane Number",                  36,  0.,        36. );
  _hPanel         = tfs->make<TH1D>( "hPanel",        "Panel Number",                   6,  0.,         6. );
  _hViewId        = tfs->make<TH1D>( "hViewId",       "View Id",                       12,  0.,        12. );
  _hStraw         = tfs->make<TH1D>( "hStraw",        "Straw Number",                  96,  0.,        96. );
  _hnBadDigis     = tfs->make<TH1D>( "hnBadDigis",    "Number Bad Digis per event",    10,   1.,       11. );
  _hnComboStraw   = tfs->make<TH1D>( "hnComboStraw",  "nCombo for Straw ComboHits",    10,   0.,       10. );
  _hnComboPanel   = tfs->make<TH1D>( "hnComboPanel",  "nCombo for Panel ComboHits",    10,   0.,       10. );
  _hnMCTracks     = tfs->make<TH1D>( "hnMCTracks",    "Number MC Tracks with Digis",   10,   0.,       10. );
  _hnMCMuons      = tfs->make<TH1D>( "hnMCMuons",     "Number MC Muons with Digis",    10,   0.,       10. );
  _hnDigisPerMuon = tfs->make<TH1D>( "hnDigsPerMuon", "Number of Digis per Muon",      50,   0.,       50. );
  _hMomentumAll   = tfs->make<TH1D>( "hMomentumAll",  "Momentum of all Muons;(MeV/c)",              100,   -20000.,   200000. );
  _hnBackground   = tfs->make<TH1D>( "hnBackground",  "Extra hits per good Muon",                    50,   0.,       50. );
  _hMomentum      = tfs->make<TH1D>( "hMomentum",     "Momentum of Muons after cuts;(MeV/c)",       100,   0.,   100000. );
  _hMomentumStart = tfs->make<TH1D>( "hMomentumStart","Start Momentum of Muons after cuts;(MeV/c)", 100,   0.,   100000. );
  _hMomentumDelta = tfs->make<TH1D>( "hMomentumDelta","Delta Momentum of Muons after cuts;(MeV/c)", 100,   0.,   100000. );
  _hnPlanes       = tfs->make<TH1D>( "hnPlanes",      "Number of planes hit by Muon",                36,   0.,       36. );
  _hWireResStraw  = tfs->make<TH1D>( "hWireResStraw", "Resolution along the wire StrawHit;(mm)",    200,   0.,      200. );
  _hWireResPanel  = tfs->make<TH1D>( "hWireResPanel", "Resolution along the wire PanelHit;(mm)",    200,   0.,      200. );
  _hTransResStraw = tfs->make<TH1D>( "hTransResStraw", "Transverse resolution StrawHit;(mm)",  60,   0.,      6. );
  _hTransResPanel = tfs->make<TH1D>( "hTransResPanel", "Transverse resolution PanelHit;(mm)",  60,   0.,      6. );

  _hUniquePanelGoodMu = tfs->make<TH1D>( "hUniquePanelGoodMu", "Unique Panel ID Good Muon",  216,  0.,       216. );
  _hPlaneGoodMu       = tfs->make<TH1D>( "hPlaneGoodMu",       "Plane Number Good Muon",      36,  0.,        36. );
  _hPanelGoodMu       = tfs->make<TH1D>( "hPanelGoodMu",       "Panel Number Good Muon",       6,  0.,         6. );
  _hViewIdGoodMu      = tfs->make<TH1D>( "hViewIdGoodMu",      "View Id Good Muon",           12,  0.,        12. );
  _hStrawGoodMu       = tfs->make<TH1D>( "hStrawGoodMu",       "Straw Number Good Muon",      96,  0.,        96. );
  _hFracPanelHits     = tfs->make<TH1D>( "hFracPanelHits",     "Frac of Straw Hits in Panel Hits, Good Muon", 101,  0., 1.01 );
  _hnMissingPanelHits = tfs->make<TH1D>( "hnMissingPanelHits", "Number Missing Straw Hits in Panel Hits, Good Muon", 20,  0., 20 );

  _hnCaloDigiAll  = tfs->make<TH1D>( "hnCaloDigiAll", "Number of Calo digis, all events",             30,   1.,       31. );
  _hnCaloDigi     = tfs->make<TH1D>( "hnCaloDigi",     "Number of Calo digis, events with good muon", 30,   1.,       31. );

  _ntWireRes = tfs->make<TNtuple>( "ntWireRes", "Wire Resolution",
                                   "istraw:ncombo:wdist:wRes:transRes");

}

bool mu2e::CosmicMuonInfo::filter(art::Event& event) {

  bool retval{false};
  auto strawDigis   = event.getValidHandle<StrawDigiCollection>( _strawDigisTag );
  auto strawDigiMCs = event.getValidHandle<StrawDigiMCCollection>( _strawDigiMCsTag );
  auto strawHits    = event.getValidHandle<StrawHitCollection>( _strawHitsTag );
  auto comboHits    = event.getValidHandle<ComboHitCollection>( _strawHitsTag );
  auto panelHits    = event.getValidHandle<ComboHitCollection>( _panelHitsTag );
  auto caloDigis    = event.getValidHandle<CaloDigiCollection>( _caloDigisTag );

  _hNStrawHits->Fill(strawHits->size());
  _hNPanelHits->Fill(panelHits->size());
  _hnCaloDigiAll->Fill(caloDigis->size());

  DigisBySim digisBySim( *strawDigiMCs );

  _hnMCTracks->Fill( digisBySim.nTracks());
  _hnMCMuons->Fill( digisBySim.nMuons());

  int sum{0};
  for ( auto const& trkinfo : digisBySim.bySim() ){
    sum += trkinfo.second.size();
  }

  if ( _diagLevel > 0 ) {
    cout << "Event: " << event.id()
         << " strawDigis:   " << strawDigis->size()
         << " strawDigisMC: " << strawDigiMCs->size()
         << " caloDigis: "    << caloDigis->size()
         << " nTracks:   "    << digisBySim.nTracks()
         << " sum:       "    << sum
         << " nGood:     "    << digisBySim.nGoodTracks(15)
         << " nBad:      "    << digisBySim.nBad()
         << endl;
  }

  _hnBadDigis->Fill( digisBySim.nBad() );

  for ( auto const& combo : *comboHits ){
    auto id0 = strawHits->at(combo.indexArray().at(0)).strawId();
    _hUniquePanel->Fill(id0.uniquePanel());
    _hPlane->Fill(id0.getPlane());
    _hPanel->Fill(id0.getPanel());
    _hViewId->Fill( viewId(id0) );
    _hStraw->Fill(id0.getStraw());
    _hWireResStraw->Fill(combo.wireRes());
    _hTransResStraw->Fill(combo.transRes());
    _hnComboStraw->Fill(combo.nCombo());
  }

  std::vector<int> strawHit2panelHit(strawHits->size(),0);

  float nt[5];
  int icombo(-1);
  for ( auto const& combo : *panelHits ){
    ++icombo;
    _hWireResPanel->Fill(combo.wireRes());
    _hTransResPanel->Fill(combo.transRes());
    _hnComboPanel->Fill(combo.nCombo());
    nt[0] = combo.strawId().getStraw();
    nt[1] = combo.nCombo();
    nt[2] = combo.wireDist();
    nt[3] = combo.wireRes();
    nt[4] = combo.transRes();
    _ntWireRes->Fill(nt);
    for ( int i : combo.indexArray() ){
      strawHit2panelHit.at(i) = icombo;
    }
  }

  for ( auto const& trkinfo : digisBySim.bySim() ){
    auto const& sim  = *trkinfo.first;
    const int nBackground = strawDigis->size() - trkinfo.second.digi_indices.size();
    const double pStart   = sim.startMomentum().vect().mag();
    const double p        = trkinfo.second.p;
    const double pDelta   = pStart-p;
    const bool isMuon = (std::abs(sim.pdgId()) == PDGCode::mu_minus);
    XYZVec momStart(sim.startMomentum().vect().x(),sim.startMomentum().vect().y(), sim.startMomentum().vect().z());
    //const double phi_start = atan2(sim.startMomentum().vect().y(),sim.startMomentum().vect().z());
    //double mag = sqrt((sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x()));
    
    //const double theta_start = acos(sim.startMomentum().vect().z()/mag);
    //const double phi_start = atan(sim.startMomentum().vect().y()/sim.startMomentum().vect().x());
    
    if ( isMuon ) {
      _hnDigisPerMuon->Fill( trkinfo.second.digi_indices.size() );
      _hMomentumAll->Fill(p);
     // _phiMC->Fill(phi_start);
      //_phiMC->SetStats(0);
      //_thetaMC->Fill(theta_start);
      //_thetaMC->SetStats(0);
     // _phiMC->SaveAs("thetaMC.root");
    }

    set<int> planes;
    for ( int i : trkinfo.second.digi_indices ){
      auto const& digi = strawDigis->at(i);
      planes.insert( digi.strawId().getPlane() );
      
    }
    if ( _diagLevel > 2 ) {
      cout << " Evt: " << event.id().event()
           << " Track: " << trkinfo.first.key()
           << " pdgid: " << sim.pdgId()
           << " hit trkinfo: " << trkinfo.second.digi_indices.size()
           << " nplanes: "    << planes.size()
           << " p: " << p
           << endl;
    }
    if ( isMuon && p  >= _cuts.pmin  &&  p  <= _cuts.pmax &&
         trkinfo.second.digi_indices.size() >= _cuts.minStrawDigis
         ){
      _hnBackground->Fill(nBackground);
      _hMomentum->Fill(p);
      _hMomentumStart->Fill(pStart);
      _hMomentumDelta->Fill(pDelta);
      _hnPlanes->Fill( planes.size() );
      //double mag = sqrt((sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x())+(sim.startMomentum().vect().x()*sim.startMomentum().vect().x()));
    
      //const double theta_cuts = acos(sim.startMomentum().vect().z()/mag);
      //const double phi_cuts = atan(sim.startMomentum().vect().y()/sim.startMomentum().vect().x());
    
      //_phiMCcuts->Fill(phi_cuts);
      //_thetaMCcuts->Fill(theta_cuts);
    

      if ( planes.size() >= _cuts.minPlanes &&
           nBackground   >= _cuts.minBackground &&
           nBackground   <= _cuts.maxBackground ){
        _hnCaloDigi->Fill(caloDigis->size());
        int nPanelHits{0};
        for ( int i : trkinfo.second.digi_indices ){
          auto id = strawDigis->at(i).strawId();
          _hUniquePanelGoodMu->Fill(id.uniquePanel());
          _hPlaneGoodMu->Fill(id.getPlane());
          _hPanelGoodMu->Fill(id.getPanel());
          _hViewIdGoodMu->Fill( viewId(id));
          _hStrawGoodMu->Fill(id.getStraw());
          if ( strawHit2panelHit.at(i) != 0 ){
            ++nPanelHits;
          }
        }
        double r = double(nPanelHits) / double(trkinfo.second.digi_indices.size());
        _hFracPanelHits->Fill(r);
        _hnMissingPanelHits->Fill( trkinfo.second.digi_indices.size() - nPanelHits);

        retval = true;
      }
    }
  }
  return retval;

}

DEFINE_ART_MODULE(mu2e::CosmicMuonInfo);
