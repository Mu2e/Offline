//
// Create printout that can be used to verify that the input and output
// data products used in event mixing are the same.  The printout
// navigates the art::Ptrs in the event and makes printout using their
// pointees. The printout does not include indices/keys because these
// may be changed by mixing.
//
// $Id: InspectDataProducts_module.cc,v 1.5 2013/10/21 20:34:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/21 20:34:14 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MixingSummary.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class InspectDataProducts : public art::EDAnalyzer {
  public:

    explicit InspectDataProducts(fhicl::ParameterSet const& pset);
    virtual ~InspectDataProducts() { }

    virtual void beginJob();
    virtual void endJob();

    // This is called for each event.
    virtual void analyze(const art::Event& e);

  private:

    // Start: run time parameters

    // Diagnostics printout level
    int _diagLevel;

    // Module label of the g4 module that made the hits.
    std::string _mixerModuleLabel;

    // Module label of the generator module that was passed as input to G4.
    std::string _generatorModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // End: run time parameters

    int _nBadG4Status;

    //
    void doStepPointMCInstance ( art::Event const&  event,
                                 std::string const& moduleLabel,
                                 std::string const& instanceName );

  };
  InspectDataProducts::InspectDataProducts(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _mixerModuleLabel(pset.get<string>("mixerModuleLabel","mixer")),
    _generatorModuleLabel(pset.get<string>("generatorModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),

    // Remaining member data
    _nBadG4Status(0){
  }

  void InspectDataProducts::beginJob(){
  }

  void InspectDataProducts::analyze(const art::Event& event) {

    art::Handle<StatusG4> g4StatusHandle;
    event.getByLabel( "g4run", g4StatusHandle);

    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel,gensHandle);
    GenParticleCollection const& gens(*gensHandle);

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(_mixerModuleLabel,simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    art::Handle<PointTrajectoryCollection> trajectoriesHandle;
    event.getByLabel(_mixerModuleLabel,trajectoriesHandle);
    PointTrajectoryCollection const& trajectories(*trajectoriesHandle);

    art::Handle<MixingSummary> mixsumHandle;
    event.getByLabel(_mixerModuleLabel,mixsumHandle);

    cout << "\n New Event: "
         << event.id()   << " | "
         << gens.size()  << " "
         << sims.size()  << " "
         << endl;


    // All StepPointMC Collections
    std::string unknown = StepInstanceName().name();
    std::vector<std::string> const& stepNames = StepInstanceName::names();
    for ( std::vector<std::string>::const_iterator i=stepNames.begin(), e=stepNames.end();
          i != e; ++i ){
      if ( *i != unknown ) {
        doStepPointMCInstance( event, _mixerModuleLabel, *i );
      }
    }

    for ( SimParticleCollection::const_iterator i=sims.begin(), e=sims.end();
          i != e; ++ i){

      SimParticleCollection::key_type key(i->first);
      SimParticle const& sim(i->second);

      if ( sim.id() != key ){
        cout << "Fubar SimParticle key match at event "
             << event.id() << " : "
             << sim.id()   << " "
             << key        << " "
             << endl;
      }

      cout << "SimParticle: "
           << sim.pdgId()            << " "
           << sim.isPrimary()        << " "
           << sim.startPosition()    << " "
           << sim.creationCode()     << " "
           << sim.stoppingCode()     << " "
           << sim.daughters().size() << " ";
      if ( sim.isPrimary() ){
        GenParticle const& gen(*sim.genParticle());
        cout << gen.time() <<  " ";
      } else {
        SimParticle const& mother(*sim.parent());
        cout << mother.startGlobalTime() << " ";
      }
      for( size_t i=0; i<sim.daughters().size(); ++i ){
        SimParticle const& dau(*sim.daughters().at(i));
        cout << dau.startGlobalTime() << " ";
      }
      cout << endl;

    } // end loop over simparticles.

    for ( GenParticleCollection::const_iterator i=gens.begin(), e=gens.end();
          i != e; ++i ){
      cout << "GenParticle: " << *i << endl;
    }

    for ( PointTrajectoryCollection::const_iterator i=trajectories.begin(), e=trajectories.end();
          i !=e ; ++i ){

      PointTrajectoryCollection::key_type key(i->first);
      PointTrajectory const&              traj(i->second);

      if ( traj.simId() != int(key.asInt()) ){
        cout << "Fubar PointTrajectory key match at event "
             << event.id()   << " : "
             << traj.simId() << " "
             << key          << " "
             << endl;
      }

      SimParticle const& sim = sims[cet::map_vector_key(traj.simId())];
      cout << " Trajectory: "
           << sim.pdgId() << " "
           << sim.startGlobalTime() << " "
           << traj.size() << " |" ;
      size_t nmax = ( traj.size() < 10 ) ? traj.size() : 10;
      for ( size_t i=0; i<nmax; ++i){
        cout << " " << traj.points().at(i).z();
      }
      cout << endl;
    }

    if ( g4StatusHandle.isValid() ){
      cout << "StatusG4: " << event.id() << " " << *g4StatusHandle << endl;
    }


    if ( mixsumHandle.isValid() ){
      MixingSummary const& mixsum(*mixsumHandle);

      std::vector<StatusG4> const& eventStatus(mixsum.eventStatus());
      art::EventIDSequence const&  eventIDs(mixsum.eventIDs());
      for ( size_t i=0; i < eventStatus.size(); ++i ){
        cout << "StatusG4: " << eventIDs.at(i) << " " << eventStatus.at(i) << endl;
      }
    }

  } // end analyze

  void InspectDataProducts::doStepPointMCInstance ( art::Event const&  event,
                                                    std::string const& moduleLabel,
                                                    std::string const& instanceName ){

    art::Handle<StepPointMCCollection> stepsHandle;
    event.getByLabel(moduleLabel,instanceName,stepsHandle);
    if ( !stepsHandle.isValid() ){
      cout << "No such stepPointMCCollection: "
           << moduleLabel <<  " "
           << instanceName << " "
           << endl;
      return;
    }

    StepPointMCCollection const& steps(*stepsHandle);
    cout << "Size: " << instanceName << " " << steps.size() << endl;

    // Loop over all hits.
    for ( size_t i=0; i<steps.size(); ++i ){

      const StepPointMC& hit = steps[i];

      SimParticle const& sim = *hit.simParticle();
      int pdgId = sim.pdgId();
      GenId genId = ( sim.fromGenerator() ) ? sim.genParticle()->generatorId() : GenId();

      cout << instanceName                     << ": "
           << hit.strawId()                    << " "
           << hit.time()                       << " "
           << hit.ionizingEdep()/CLHEP::keV    << " "
           << hit.stepLength()                 << " | "
           << pdgId                            << " "
           << hit.endProcessCode()             << " "
           << sim.startGlobalTime()            << " "
           << sim.startMomentum().vect().mag() << " "
           << genId
           << endl;

    } // end loop over hits.

  } // end doStepPointMCInstance

  void InspectDataProducts::endJob(){
    cout << "InspectDataProducts::endJob Number of events skipped "
         << "due to G4 completion status: "
         << _nBadG4Status
         << endl;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::InspectDataProducts);
