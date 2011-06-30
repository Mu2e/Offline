//
// Read the mixed events.
//
// $Id: MixAnalyzer_module.cc,v 1.1 2011/06/30 04:38:47 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/30 04:38:47 $
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

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

namespace mu2e {

  class MixAnalyzer : public art::EDAnalyzer {
  public:
    explicit MixAnalyzer(fhicl::ParameterSet const& pset){}
    virtual ~MixAnalyzer() { }

    void analyze( art::Event const& e);

    void endJob();

  private:

  };

  void
  MixAnalyzer::analyze(art::Event const& event) {

    //mf::LogInfo log("Mix");

    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel("mixer",gensHandle);
    GenParticleCollection const& gens = *gensHandle;

    art::Handle<MixingSummary> sumHandle;
    event.getByLabel("mixer",sumHandle);
    MixingSummary const& summary = *sumHandle;

    for ( size_t i=0; i != summary.eventIDs().size(); ++i ){
      cerr << "Mixed event: " << summary.eventIDs()[i] << endl;
    }

    for ( size_t i=0; i!=gens.size(); ++i){
      cerr << "Analyze: gen "
           << i << " "
           << gens[i]
           << endl;
    }

  } // end of ::analyze.

  void MixAnalyzer::endJob(){
    cerr << "All done with analyzer: " << endl;
  }
}

using mu2e::MixAnalyzer;
DEFINE_ART_MODULE(MixAnalyzer);
