//
// Illustrate a bug in the mixing infrastructure:
//
// $Id: MixBug01_module.cc,v 1.2 2011/10/28 18:47:06 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:06 $
//
// Contact person Rob Kutschke.
//
// See comments in MixMCEvents_module.cc for details of how this module works.

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"

// Includes from art
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "canvas/Utilities/InputTag.h"

// Other third party includes
#include "boost/noncopyable.hpp"

// C++ includes
#include <vector>

namespace mu2e {
  class MixBug01Detail;
  typedef art::MixFilter<MixBug01Detail> MixBug01;
}

class mu2e::MixBug01Detail : private boost::noncopyable {

public:

  MixBug01Detail(fhicl::ParameterSet const &p,
                 art::MixHelper &helper);

  void startEvent();

  size_t nSecondaries();

  void processEventIDs(art::EventIDSequence const &seq);

  void finalizeEvent(art::Event &t);

  bool 
  aggregateGens( std::vector<mu2e::GenParticleCollection const *> const& in,
                 mu2e::GenParticleCollection&                            out,
                 art::PtrRemapper const & );

private:

  // Count number of events.
  int count_;

  // Number of secondaries for this event.
  size_t actual_;

};

mu2e::MixBug01Detail::
MixBug01Detail(fhicl::ParameterSet const &pSet,
               art::MixHelper &helper):
  count_(-1){

  helper.declareMixOp
    ( art::InputTag("generate",""),
      &MixBug01Detail::aggregateGens, *this );

}

void
mu2e::MixBug01Detail::
startEvent() {
  std::cerr << "\nStart event: " << ++count_ << std::endl;
}

size_t
mu2e::MixBug01Detail::
nSecondaries() {

  // Step on the landmine in event 2.
  int n = (count_ == 2) ? 0 : 2;

  actual_ = n;
  std::cerr << "nSecondaries: " << actual_ << std::endl;
  return n;
}

void
mu2e::MixBug01Detail::
processEventIDs(art::EventIDSequence const &seq) {

  std::string tag = ( actual_ == seq.size() ) ? "" : "fubar";
  std::cerr << "processEventIDs: " 
            << actual_ << " "
            << seq.size() << " "
            << tag
            << std::endl;

  for ( size_t i=0; i<seq.size(); ++i){
    std::cerr << "   Ids: "
              << i     << " : "
              << seq.at(i) << " "
              << "\n";
  }
}

void
mu2e::MixBug01Detail::
finalizeEvent(art::Event &e) {
  std::cerr << "Finalize ... " << e.event() << std::endl;
}

bool
mu2e::MixBug01Detail::
aggregateGens( std::vector< mu2e::GenParticleCollection const *> const& in,
               mu2e::GenParticleCollection&                             out,
               art::PtrRemapper const & ){

  std::string tag = ( actual_ == in.size() ) ? "" : "fubar";
  std::cerr << "aggregateGens: " 
            << actual_ << " "
            << in.size() << " "
            << tag
            << std::endl;

  for ( std::vector< mu2e::GenParticleCollection const *>::const_iterator i=in.begin(), e=in.end();
        i != e ; ++i ){
    mu2e::GenParticleCollection const& gens(**i);
    std::cerr << "  Collection size: " << gens.size() << std::endl;
    for ( mu2e::GenParticleCollection::const_iterator g=gens.begin(), ge=gens.end();
          g !=ge; ++g){
      std::cerr << "      " << *g << std::endl;
    }
            
  }

  return true;
} // end mu2e::MixBug01Detail::aggregateGens


DEFINE_ART_MODULE(mu2e::MixBug01);
