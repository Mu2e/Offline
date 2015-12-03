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
    for ( auto const& dev : tracker.getPlanes() ){
      for ( auto const& sec : dev.getPanels() ){
        StrawId sid( sec.id(), 0, 0 );
        Straw const& straw = sec.getStraw(sid);
        double phi  = straw.direction().phi();
        double z    = straw.getMidPoint().z() - dev.origin().z();
        double phi1 = phi/M_PI*180.;
        cout << "panel: "
             << sec.id()      << " "
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
