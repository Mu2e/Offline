/*
 *
 * Two different ways of looking at provenances:
 *
 * 1) Get a handle to a named container and then look at its provenance.
 *    This also shows how to access the ParameterSet used to create the 
 *    Data Product.
 * 2) Get a list of all handles in the event and look at each of them.

 * $Id: Ex03InspectProvenance_plugin.cc,v 1.3 2011/05/17 15:36:00 greenc Exp $
 * $Author: greenc $
 * $Date: 2011/05/17 15:36:00 $
 *  
 * Original author Rob Kutschke
 *
 */

// C++ includes
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

// Mu2e includes
#include "ToyDP/inc/ToyHitCollection.hh"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class Ex03InspectProvenance : public art::EDAnalyzer {
  public:
    explicit Ex03InspectProvenance(fhicl::ParameterSet const& pset):
      _done(false){ 
    }
    virtual ~Ex03InspectProvenance() { }

    void analyze(const art::Event& e, art::EventSetup const&);

  private:

    // Only do this once.
    bool _done;

  };
  
  void
  Ex03InspectProvenance::analyze(const art::Event& evt, art::EventSetup const&) {

    // Only do this once.
    if ( _done ) return;

    // Name of the module that created the hits of interest;
    static const string creatorName("ex02hitmaker");

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<ToyHitCollection> handle;
    evt.getByLabel(creatorName,handle);
    
    // Only say we are done after we actually found something in the event.
    _done = true;
    
    // Get the provenance of the data product.
    const art::Provenance& prov = *(handle.provenance());

    // Print the provenance and keep the message open for additions.
    mf::LogInfo log("ProvenanceInfo");
    log << "Provenance of the hits in this event: \n"
         << prov 
         << "\n";
    
    // Extract the parameter set IDs from the provenance.
    const set<art::ParameterSetID>&  ids = prov.psetIDs();
    log << "Number of Parameter set IDs: " << ids.size() << "\n";

    // Parameter set Registry singleton.
    art::pset::Registry* reg = art::pset::Registry::instance();

    // Loop over all parameter set IDs.
    set<art::ParameterSetID>::const_iterator b = ids.begin();
    set<art::ParameterSetID>::const_iterator e = ids.end();
    for ( ; b!=e; ++b){
      fhicl::ParameterSet result;

      if (!reg->getMapped(*b, result)){
        log << "Could not get the parameter set... \n";
      }else{
        log << result;

        // We know that there is only one provenance and that it has
        // a parameter named minPulseHeight.
        double minPulseHeight = result.get<double>("minPulseHeight");
        log << "\nminimum pulse height: " << minPulseHeight << "\n";
      }
    }

    // Get all of the provenances in this event.
    vector<art::Provenance const*> Prov;
    evt.getAllProvenance( Prov );
    log << "\nNumber of provenances in this event: " << Prov.size();

    //Loop over all of the provenances in this event.
    for ( vector<art::Provenance const *>::size_type i=0;
          i<Prov.size(); ++i ){
      const art::Provenance&  prov = *Prov[i];
      log  << "\nNext Provenance: \n"
           << prov; 
    }
  }
}

using mu2e::Ex03InspectProvenance;
DEFINE_ART_MODULE(Ex03InspectProvenance);
