//
// A stripped down version of ReadStrawHit that only looks at reco info
//
// Original author Rob Kutschke.
//

// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1F.h"
#include "TNtuple.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class ReadStrawHitReco : public art::EDAnalyzer {
  public:
    explicit ReadStrawHitReco(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;
    virtual void analyze( art::Event const& e) override;

  private:

    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Input tag describing which hits to look at
    art::InputTag _hitsTag;

    // Some diagnostic histograms.
    TH1F* _hStrawNumber    = nullptr;
    TH1F* _hLayerNumber    = nullptr;
    TH1F* _hPanelNumber    = nullptr;
    TH1F* _hPlaneNumber    = nullptr;
    TH1F* _hHitTime        = nullptr;
    TH1F* _hHitDeltaTime   = nullptr;
    TH1F* _hHitEnergy      = nullptr;
    TH1F* _hNHits          = nullptr;
    TH1F* _hNHits1         = nullptr;
    TH1F* _hNHitsPerWire   = nullptr;
    TH1F* _hFractionOnWire = nullptr;
    TH1F* _hTDivOK         = nullptr;
    TNtuple* _ntup         = nullptr;

    ProditionsHandle<StrawResponse> _strawResponse_h;
  };

} // end of namespace mu2e

mu2e::ReadStrawHitReco::ReadStrawHitReco(fhicl::ParameterSet const& pset):
  art::EDAnalyzer(pset),
  _diagLevel(pset.get<int>("diagLevel",0)),
  _maxFullPrint(pset.get<int>("maxFullPrint",0)),
  _hitsTag(pset.get<std::string>("hitsTag")){
}

void mu2e::ReadStrawHitReco::beginJob(){

  if ( _diagLevel > 0 ) {
    cout << "ReadStrawHitReco Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint
         << endl;
  }

  art::ServiceHandle<art::TFileService> tfs;

  _hStrawNumber    = tfs->make<TH1F>( "hStrawNumber",    "Straw Number in Panel",          100,   0.,   100. );
  _hLayerNumber    = tfs->make<TH1F>( "hLayerNumber",    "Straw Number in Panel",            4,   0.,     2. );
  _hPanelNumber    = tfs->make<TH1F>( "hPanelNumber",    "Straw Number in Panel",           12,   0.,     6. );
  _hPlaneNumber    = tfs->make<TH1F>( "hPlaneNumber",    "Straw Number in Panel",           50,   0.,    50. );
  _hHitTime        = tfs->make<TH1F>( "hHitTime",        "Hit Time (ns)",                  200,   0.,  2000. );
  _hHitDeltaTime   = tfs->make<TH1F>( "hHitDeltaTime",   "Hit Delta Time (ns)",             80, -20.0,   20. );
  _hHitEnergy      = tfs->make<TH1F>( "hHitEnergy",      "Hit Energy (keV)",                75,   0.,    15. );
  _hNHits          = tfs->make<TH1F>( "hNHits",          "Number of straw hits",           500,   0.,   500. );
  _hNHits1         = tfs->make<TH1F>( "hNHits1",         "Number of straw hits",           200,   0., 12000. );
  _hNHitsPerWire   = tfs->make<TH1F>( "hNHitsPerWire",   "Number of hits per straw",        10,   0.,    10. );
  _hFractionOnWire = tfs->make<TH1F>( "hFractionOnWire", "Fractional Position along wire", 100,  -1.5,    1.5 );
  _hTDivOK         = tfs->make<TH1F>( "hTDivOK",         "Time division OK",                 4,   0.,     2.  );
  _ntup            = tfs->make<TNtuple>( "ntup", "Straw Hit ntuple",
                                         "evt:index:sid:lay:did:panel:time:dtime:eDep:hitx:hity:hitz:frac:tdivOK");
  }

void mu2e::ReadStrawHitReco::analyze(art::Event const& evt) {

  static int ncalls(0);
  ++ncalls;

  const Tracker& tracker = *GeomHandle<Tracker>();

  auto const& hits = *evt.getValidHandle<StrawHitCollection>( _hitsTag );
  if ( _diagLevel > 1 ) {
    cout << "ReadStrawHitReco: Total number of straw hits = " << hits.size() << endl;
  }
  _hNHits->Fill(hits.size());
  _hNHits1->Fill(hits.size());

  auto const& srep = _strawResponse_h.get(evt.id());

  // Counter for number of hits on each wire.
  std::map<StrawId,int> nhperwire;

  for ( StrawHit const& hit : hits ) {

    _hHitTime ->Fill(hit.time());
    _hHitDeltaTime->Fill(hit.dt());
    _hHitEnergy->Fill(hit.energyDep()*1000.0);

    StrawId sid = hit.strawId();

    unsigned const id      = sid.asUint16();
    Straw const& str       = tracker.getStraw(sid);
    unsigned const lno     = sid.getLayer();
    unsigned const plno    = sid.getPlane();
    unsigned const pnno     = sid.getPanel();

    _hStrawNumber->Fill(sid.getStraw());
    _hLayerNumber->Fill(lno);
    _hPanelNumber->Fill(pnno);
    _hPlaneNumber->Fill(plno);

    double halfLen                  = str.halfLength();

    // Calculate the hit position; it's a point on the wire
    // (for now wire is modeled as a straight line).
    double dw, dwerr, halfpv;
    bool td = srep.wireDistance(str,hit.energyDep(),hit.dt(),dw,dwerr,halfpv);
    _hTDivOK->Fill( td );

    // Fractional length along the wire.
    // Only fill the histogram if time division info is useful.
    double frac = dw/halfLen;
    CLHEP::Hep3Vector pos = str.getMidPoint()+dw*str.getDirection();
    if ( td ) {
      _hFractionOnWire->Fill(frac);
    }

    // Fill the ntuple:
    float nt[_ntup->GetNvar()];
    nt[0]  = evt.id().event();
    nt[1]  = id;
    nt[2]  = sid.getStraw();
    nt[3]  = lno;
    nt[4]  = plno;
    nt[5]  = pnno;
    nt[6]  = hit.time();
    nt[7]  = hit.dt();
    nt[8]  = hit.energyDep();
    nt[9]  = pos.getX();
    nt[10] = pos.getY();
    nt[11] = pos.getZ();
    nt[12] = frac;
    nt[13] = td;

    _ntup->Fill(nt);

    // Calculate number of hits per wire
    ++nhperwire[hit.strawId()];

    if ( (int(evt.id().event()) < _maxFullPrint) && (_diagLevel > 3) ) {
      cout << "ReadStrawHitReco: "
           << evt.id().event()      << " #"
           << sid                   << " "
           << hit.time()            << " "
           << hit.dt()              << " "
           << hit.energyDep()
           << endl;
    }

  } // end loop over hits.

  for( std::map<StrawId,int>::iterator it=nhperwire.begin(); it!= nhperwire.end(); ++it ) {
    _hNHitsPerWire->Fill(it->second);
  }

} // end of ::analyze.

DEFINE_ART_MODULE(mu2e::ReadStrawHitReco);
