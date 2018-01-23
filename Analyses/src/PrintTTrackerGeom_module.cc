//
// Print the information about the TTracker
//
// $Id: PrintTTrackerGeom_module.cc,v 1.1 2013/12/20 20:05:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/12/20 20:05:12 $
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

  };

  PrintTTrackerGeom::PrintTTrackerGeom(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset){
  }

  void PrintTTrackerGeom::analyze(const art::Event& ){}

  void PrintTTrackerGeom::beginRun(const art::Run& run){

    TTracker const& tracker(*GeomHandle<TTracker>());

    cout << "Tracker: " << tracker.nPlanes() << endl;
    // for ( auto const& pln : tracker.getPlanes() ){
    for ( size_t i=0; i!= tracker.nPlanes(); ++i){
      const auto& pln = tracker.getPlane(i);
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

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PrintTTrackerGeom);
