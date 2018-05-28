//
// Print the information about the TTracker
//
// Original author Rob Kutschke
//

#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class PrintTTrackerGeom : public art::EDAnalyzer {
  public:

    explicit PrintTTrackerGeom(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e) override;

    void beginRun ( const art::Run& r) override;

  private:
    int _diagLevel;

  };

  PrintTTrackerGeom::PrintTTrackerGeom(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset),
    _diagLevel(pset.get<int>("diagLevel",0)){
  }

  void PrintTTrackerGeom::analyze(const art::Event& ){}

  void PrintTTrackerGeom::beginRun(const art::Run& run){

    TTracker const& tracker(*GeomHandle<TTracker>());

    cout << "Tracker: " << tracker.nPlanes() << endl;
    for ( auto const& pln : tracker.getPlanes() ){
      for ( auto const& pnl : pln.getPanels() ){
        StrawId sid( pnl.id() ); // first straw id is equal to its panel id
        Straw const& straw = pnl.getStraw(sid);
        double phi  = straw.direction().phi();
        double z    = straw.getMidPoint().z() - pln.origin().z();
        double phi1 = phi/M_PI*180.;
        cout << "panel: "
             << pnl.id()      << " "
             << sid           << " "
             << straw.index() << " : "
             << z             << " "
             << phi1
             << endl;
      }
    }

    if ( _diagLevel == 0 ) return;

    auto const& straws = tracker.getAllStraws();
    cout << "Straws: " << straws.size() << "  " << StrawId::_end << endl;
    size_t n(0);
    for ( Straw const& straw : straws ){

      cout << "Straw1: "
           << n++ << "  "
           << straw.index().asInt() << " "
           << straw.id().asUint16() << " "
           << straw.id() << "  |  "
           << straw.getMidPoint() << " "
           << straw.direction()
           << endl;
    }

    n=0;
    for ( uint16_t ipln=0; ipln < StrawId::_nplanes ; ++ipln ){
      for ( uint16_t ipan=0; ipan < StrawId::_npanels ; ++ipan ){
        for ( uint16_t istr=0; istr < StrawId::_nstraws; ++istr ){
          StrawId sid( ipln, ipan, istr);
          Straw const& straw = tracker.getStraw( sid );
          cout << "Straw2: "
               << n++ << " "
               << sid.asUint16() << " "
               << sid  << "  | ";
          for ( auto nid : straw.nearestNeighboursById() ){
            cout << " " << nid;
          }
          for ( auto idx : straw.nearestNeighboursByIndex() ){
            cout << " " << tracker.getStraw(idx).id();
          }
          cout << " "
               << straw.getMidPoint() << " "
               << straw.direction()
               << endl;
        }
      }
    }
    cout << "Counted Straws: " << n << endl;

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PrintTTrackerGeom);
