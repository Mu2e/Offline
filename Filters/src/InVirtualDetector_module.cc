//
// Given the id number of a virtual detector and a list of pdgIds,
// select events in which at least on of the particles has a step in
// that virtual detector.
//
// $Id: InVirtualDetector_module.cc,v 1.2 2013/05/30 18:40:35 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/30 18:40:35 $
//
// Contact person Rob Kutschke.
//

// Mu2e includes.
#include "Mu2eUtilities/inc/DiagnosticsG4.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"


// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// Root includes
#include "TH1F.h"
#include "TNtuple.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

namespace mu2e {

  class InVirtualDetector : public art::EDFilter {

  public:
    explicit InVirtualDetector(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor

    bool filter( art::Event& event);

    virtual void beginJob();
    virtual void endJob();

  private:

    // Module label of the module that made the StepPointMCCollections.
    std::string g4ModuleLabel_;

    // Instance names of the StepPointMCCollections.
    StepInstanceName instance_;

    // Id of the virtual detector we care about.
    VirtualDetectorId vdId_;

    // Only count particles on this list of pdgIds.
    std::vector<int> pdgIds_;

    // Control printout
    int verbosity_;
    int maxPrint_;
    int nEvents_;

    // Histograms
    DiagnosticsG4 diagnostics_;
    TH1F*         hnSteps_;
    TH1F*         hnGood_;

    // Number of events that pass the filter.
    int nPassed_;

    // Helper functions.
    void printConfig( ostream& );

  };

  InVirtualDetector::InVirtualDetector(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    g4ModuleLabel_(pset.get<string>("g4ModuleLabel")),
    instance_(StepInstanceName::virtualdetector),
    vdId_(pset.get<int>("vdId")),
    pdgIds_(pset.get<vector<int> >("pdgIds")),
    verbosity_(pset.get<int>("verbosity",0)),
    maxPrint_(pset.get<int>("maxPrint",0)),
    nEvents_(-1),
    hnSteps_(0),
    hnGood_(0),
    nPassed_(0){

    if ( verbosity_ > 0 ) printConfig(cout);
  }

  void InVirtualDetector::beginJob(){

    diagnostics_.book("BadG4Status");

    art::ServiceHandle<art::TFileService> tfs;
    hnSteps_ = tfs->make<TH1F>( "nSteps", "Number of steps in virtual detectors",   200, 0., 200.  );
    hnGood_  = tfs->make<TH1F>( "nGood",  "Number of steps passing filter criteria", 20, 0.,  20.  );

  }

  bool
  InVirtualDetector::filter(art::Event& event) {

    ++nEvents_;

    // Accept only events with good status from G4.
    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( g4ModuleLabel_, g4StatusHandle);
    StatusG4 const& g4Status = *g4StatusHandle;

    if ( g4Status.status() > 1 ) {

      // Diagnostics for rejected events.
      diagnostics_.fillStatus( g4Status);
      return false;
    }


    // Get the steps and document their number.
    art::Handle<StepPointMCCollection> stepsHandle;
    event.getByLabel(g4ModuleLabel_, instance_.name(), stepsHandle);
    StepPointMCCollection const& steps(*stepsHandle);

    hnSteps_->Fill(steps.size());
    if ( verbosity_ > 1 ) {
      cout << "InVirtualDetector::filter Number of steps: "
           << steps.size()
           << endl;
    }

    // Number of steps that pass the filter criteria
    int nGood(0);

    for ( StepPointMCCollection::const_iterator i=steps.begin();
          i != steps.end(); ++i ){

      StepPointMC const& step(*i);
      SimParticle const& sim(*step.simParticle());

      // Apply filter criteria
      bool ok(false);
      if ( step.volumeId() == vdId_ ){
        if ( find( pdgIds_.begin(), pdgIds_.end(), sim.pdgId() )
             != pdgIds_.end() ){
          ok = true;
          ++nGood;
        }
      }

      if ( verbosity_ > 2 && nEvents_ < maxPrint_ ){
        cout << "InVirtualDetector::filter step: "
             << i-steps.begin() << " "
             << step.volumeId() << " "
             << sim.pdgId() << " "
             << ok << " "
             << nGood
             << endl;
      }
    }

    // Document the result.
    hnGood_->Fill(nGood);
    if ( verbosity_ >  1 && nEvents_ < maxPrint_ ) {
      cout << "InVirtualDetector::filter number of accepted steps: "
           << nGood
           << endl;
    }

    if ( nGood > 0 ){
      nPassed_++;
      return true;
    }

    return false;

  } // end of ::filter

  void InVirtualDetector::endJob() {
    cout << "InVirtualDetector::endJob: Number of events passing the filter: "
         << nPassed_
         << "\n";
  }

  void InVirtualDetector::printConfig( ostream& out ){

    out << "\nInVirtualDetector configuration: \n"
        << "   Virtual detector to select: "
        << vdId_
        << endl;

    GlobalConstantsHandle<ParticleDataTable> pdt;

    out << "Particles to select: ";

    for ( size_t i=0; i<pdgIds_.size(); ++i ){

      int pdgId = pdgIds_[i];

      string space = (i == 0) ? " " : ", ";

      ParticleDataTable::maybe_ref particle = pdt->particle(pdgId);
      if ( ! particle ) {
        throw cet::exception("CONFIG")
          << "Do not recognize this pdg Id number: " << pdgId << "\n";
      }
      out << space << pdgIds_[i] << " " << particle.ref().name();

    }
    out << " ]" << endl;

  } // end InVirtualDetector::printConfig

}

using mu2e::InVirtualDetector;
DEFINE_ART_MODULE(InVirtualDetector);
