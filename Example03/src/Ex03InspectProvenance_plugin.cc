/*
 *
 * Two different ways of looking at provenances:
 *
 * 1) Get a handle to a named container and then look at its provenance.
 *    This also shows how to access the ParameterSet used to create the 
 *    Data Product.
 * 2) Get a list of all handles in the event and look at each of them.

 * $Id: Ex03InspectProvenance_plugin.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
 * $Author: kutschke $
 * $Date: 2009/09/30 22:57:47 $
 *  
 * Original author Rob Kutschke
 *
 */

// C++ includes
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

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
  class Ex03InspectProvenance : public edm::EDAnalyzer {
  public:
    explicit Ex03InspectProvenance(edm::ParameterSet const& pset):
      _done(false){ 
    }
    virtual ~Ex03InspectProvenance() { }

    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Only do this once.
    bool _done;

  };
  
  void
  Ex03InspectProvenance::analyze(const edm::Event& evt, edm::EventSetup const&) {

    // Only do this once.
    if ( _done ) return;

    // Name of the module that created the hits of interest;
    static const string creatorName("ex02hitmaker");

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<ToyHitCollection> handle;
    evt.getByLabel(creatorName,handle);
    
    // Only say we are done after we actually found something in the event.
    _done = true;
    
    // Get the provenance of the data product.
    const edm::Provenance& prov = *(handle.provenance());

    // Print the provenance and keep the message open for additions.
    edm::LogInfo log("ProvenanceInfo");
    log << "Provenance of the hits in this event: \n"
	 << prov 
	 << "\n";
    
    // Extract the parameter set IDs from the provenance.
    const set<edm::ParameterSetID>&  ids = prov.psetIDs();
    log << "Number of Parameter set IDs: " << ids.size() << "\n";

    // Parameter set Registry singleton.
    edm::pset::Registry* reg = edm::pset::Registry::instance();

    // Loop over all parameter set IDs.
    set<edm::ParameterSetID>::const_iterator b = ids.begin();
    set<edm::ParameterSetID>::const_iterator e = ids.end();
    for ( ; b!=e; ++b){
      edm::ParameterSet result;

      if (!reg->getMapped(*b, result)){
	log << "Could not get the parameter set... \n";
      }else{
	log << result;

	// We know that there is only one provenance and that it has
	// a parameter named minPulseHeight.
	double minPulseHeight = result.getParameter<double>("minPulseHeight");
	log << "\nminimum pulse height: " << minPulseHeight << "\n";
      }
    }

    // Get all of the provenances in this event.
    vector<edm::Provenance const*> Prov;
    evt.getAllProvenance( Prov );
    log << "\nNumber of provenances in this event: " << Prov.size();

    //Loop over all of the provenances in this event.
    for ( vector<edm::Provenance const *>::size_type i=0;
	  i<Prov.size(); ++i ){
      const edm::Provenance&  prov = *Prov[i];
      log  << "\nNext Provenance: \n"
	   << prov; 
    }

  }
}


using mu2e::Ex03InspectProvenance;
DEFINE_FWK_MODULE(Ex03InspectProvenance);
