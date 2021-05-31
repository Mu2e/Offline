//
// Check self consistency of all SimParticleCollections in the event.
//
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <set>
#include <cmath>
#include <algorithm>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "MCDataProducts/inc/SimParticleCollection.hh"

// Root includes.
#include "TH1F.h"

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  class SimParticleCheck00 : public art::EDAnalyzer {
  public:

    explicit SimParticleCheck00(fhicl::ParameterSet const& pset);

    virtual void beginJob( );
    virtual void endJob  ( );

    virtual void analyze ( const art::Event& event);


  private:

    int nSims;
    int nBornBeforeParent;
    set<ProcessCode> badCodes;

    TH1F* _hNCollections;
    TH1F* _hDeltaTAll;
    TH1F* _hDeltaTBad;

  };

  SimParticleCheck00::SimParticleCheck00(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    nSims(0),
    nBornBeforeParent(0),
    badCodes(),

    // Histograms
    _hNCollections(0),
    _hDeltaTAll(0),
    _hDeltaTBad(0)
  {}

  void SimParticleCheck00::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _hNCollections = tfs->make<TH1F>( "hNCollections", "Number of SimParticle Collections",    20,      0.,    20. );
    _hDeltaTAll    = tfs->make<TH1F>( "hDeltaTAll",    "Delta Time (Particle-Parent); [ns]",  200,  -3000.,  3000. );
    _hDeltaTBad    = tfs->make<TH1F>( "hDeltaTBad",    "Delta Time (Particle-Parent); [ns]",  200,  -3000.,     0. );

  } // end beginJob

  void SimParticleCheck00::analyze(const art::Event& event) {
    typedef std::vector< art::Handle<SimParticleCollection> > HandleVector;
    HandleVector allSims = event.getMany<SimParticleCollection>();

    _hNCollections->Fill(allSims.size());

    // See if any children are created before the parent!
    for ( HandleVector::const_iterator i=allSims.begin(), e=allSims.end(); i != e; ++i ){
      const SimParticleCollection& sims(**i);
      nSims += sims.size();
      for ( SimParticleCollection::const_iterator j=sims.begin(), je=sims.end(); j != je; ++j ){
        const SimParticle& sim = j->second;
        if ( sim.isSecondary() ){
          double t  = sim.startGlobalTime();
          double t0 = sim.parent()->startGlobalTime();
          double dt = t-t0;
          _hDeltaTAll->Fill(dt);
          if ( t < t0 ){
            ++nBornBeforeParent;
            badCodes.insert(sim.creationCode());
            _hDeltaTBad->Fill(dt);
          }
        }
      }
    }

  } // end analyze

  void SimParticleCheck00::endJob(){

    cout << "\nNumber of SimParticles in this file:                      " << nSims             << endl;
    cout << "Number of particles born before their parents in this file: " << nBornBeforeParent << endl;
    if ( nBornBeforeParent == 0 ) return;

    double ratio = ( nSims > 0 ) ? double(nBornBeforeParent)/double(nSims) : 0.;
    cout << "Ratio:                                                      " << ratio << endl;

    cout << "\nNumber of bad process codes that give daughters born before parent: " << badCodes.size() << endl;

    for ( set<ProcessCode>::const_iterator i=badCodes.begin(), e=badCodes.end(); i!=e; ++i){
      cout << "ProcessCode for SimParticle born before parent: " << *i << endl;
    }
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleCheck00);
