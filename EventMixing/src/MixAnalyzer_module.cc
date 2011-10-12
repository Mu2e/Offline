//
// Read the mixed events.
//
// $Id: MixAnalyzer_module.cc,v 1.2 2011/10/12 20:12:09 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/10/12 20:12:09 $
//
// Contact person Rob Kutschke.
//

// Mu2e includes.
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MixingSummary.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/Selector.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

namespace mu2e {

  class MixAnalyzer : public art::EDAnalyzer {
  public:
    explicit MixAnalyzer(fhicl::ParameterSet const& pset){}
    virtual ~MixAnalyzer() { }

    void analyze( art::Event const& e);

  private:

  };

  void
  MixAnalyzer::analyze(art::Event const& event) {

    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector("");

    typedef std::vector< art::Handle<GenParticleCollection> > gensHandleVector;

    // Get all of the tracker StepPointMC collections from the event:
    gensHandleVector gensHandles;
    event.getMany( selector, gensHandles);

    cerr << "\nGenerated particles: " << endl;
    for ( gensHandleVector::const_iterator i=gensHandles.begin(), e=gensHandles.end();
          i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      cerr << "   " << prov.branchName() << endl;
      GenParticleCollection const& gens(**i);
      for ( size_t i=0; i!=gens.size(); ++i){
        cerr << "    "
             << i  << " "
             << gens[i]
             << endl;
      }
    }

    typedef std::vector< art::Handle<MixingSummary> > sumsHandleVector;

    sumsHandleVector sumsHandles;
    event.getMany( selector, sumsHandles);
    cerr << "\nMixing Summary: " << endl;
    for ( sumsHandleVector::const_iterator i=sumsHandles.begin(), e=sumsHandles.end();
          i != e; ++i ){
      art::Provenance const& prov(*(i->provenance()));
      cerr << "Branchname: " << prov.branchName() << endl;
      MixingSummary const& summary(**i);
      cerr << summary << endl;
    }

  } // end of ::analyze.

}

using mu2e::MixAnalyzer;
DEFINE_ART_MODULE(MixAnalyzer);
