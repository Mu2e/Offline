//
// Inspect the output of the StoppedParticleReactionGunN module.
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "MCDataProducts/inc/GenParticleCollections.hh"

using namespace std;

namespace mu2e {

  class PrintReactionN : public art::EDAnalyzer {
  public:

    explicit PrintReactionN(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& e);

  private:

    art::InputTag gensTag_;

    // Copy of the most recent stash data product.
    GenParticleCollections stash_;

    // The finger points to the next free spot in the stash.
    size_t finger_     = 0;
    bool   firstEvent_ = true;

  };

}  // end namespace mu2e

namespace {

  bool bitForBitEqual( mu2e::GenParticle const& l, mu2e::GenParticle const& r){
    if ( l.pdgId()       != r.pdgId()       ) return false;
    if ( l.generatorId() != r.generatorId() ) return false;
    if ( l.time()        != r.time()        ) return false;
    if ( l.properTime()  != r.properTime()  ) return false;
    if ( l.position()    != r.position()    ) return false;
    if ( l.momentum()    != r.momentum()    ) return false;
    return true;
  }
}


mu2e::PrintReactionN::PrintReactionN(fhicl::ParameterSet const& pset) :
  art::EDAnalyzer(pset)
  , gensTag_( pset.get<std::string>("gensTag") ){
}

void mu2e::PrintReactionN::analyze(const art::Event& event) {

  auto one = event.getValidHandle<GenParticleCollection>( gensTag_ );
  auto all = event.getValidHandle<GenParticleCollections>( gensTag_ );

  // On the first event, the stash must not be empty.
  if ( firstEvent_ ){
    firstEvent_ = false;
    if ( all->size() == 0 ){
      throw cet::exception("GENERATOR")
        << "Illegal event: empty GenParticleCollections on first event !\n";
    }
  }

  // Every every event must have a non-empty standard GenParticleCollection
  if ( one->size() == 0 ){
      throw cet::exception("GENERATOR")
        << "Illegal event: empty GenParticleCollection!\n";
  }

  // New stash is in the event; save a copy and reset the finger to the front of the stash.
  if ( all->size() > 0 ){
    finger_=0;
    stash_ = *all;
  }

  // Find entry in the stash that is supposed to match this event.
  size_t const oldFinger = finger_;
  GenParticleCollection const& match = stash_.at(finger_++);

  bool isEqual = bitForBitEqual( one->at(0),  match.at(0) );

  cout << "Event "
       << event.id().event() << " | "
       << one->size()   << " "
       << all->size()   << " "
       << stash_.size() << " "
       << oldFinger     << " "
       << match.size();
  if ( isEqual ){
    cout << " Exact match" << endl;
  } else{
    cout << " Do not match" << endl;
    cout << "      " << one->at(0) << endl;
    cout << "      " << match.at(0) << endl;
  }

}

DEFINE_ART_MODULE(mu2e::PrintReactionN);
