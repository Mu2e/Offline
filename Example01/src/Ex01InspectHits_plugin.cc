/*------------------------------------------------------

  Look at some overly simplified hits that are in the event.

  $Id: Ex01InspectHits_plugin.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
  $Author: kutschke $
  $Date: 2009/09/30 22:57:47 $
   
  Original author Rob Kutschke

 --------------------------------------------------------*/

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

// Mu2e includes.
#include "ToyDP/inc/ToyHitCollection.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class Ex01InspectHits : public edm::EDAnalyzer {
  public:
    explicit Ex01InspectHits(edm::ParameterSet const& pset) : 
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _nAnalyzed(0),
      _hist1(0),
      _hist2(0),
      _hist3(0),
      _hist4(0),
      _messageCategory("ToyHitInfo"){
    }
    virtual ~Ex01InspectHits() { }

    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    virtual void beginRun(edm::Run const &r, 
			  edm::EventSetup const& eSetup );

    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lblock, 
				      edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms to be filled.
    TH1F* _hist1;
    TH2F* _hist2;
    TH1F* _hist3;
    TH1F* _hist4;

    // A category for the error logger.
    const std::string _messageCategory;

  };

  void Ex01InspectHits::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    // Create some 1D histograms.
    _hist1 = tfs->make<TH1F>( "hist1", "Radius of Hits", 100,  0., 50. );
    _hist3 = tfs->make<TH1F>( "hist3", "Pulse Height", 100,  0., 20. );
    _hist4 = tfs->make<TH1F>( "hist4", "Hits per Event", 10,  0., 10. );

    // Create a subdirectory and create a 2D histogram in the subdirectory.
    edm::TFileDirectory tfdir = tfs->mkdir( "subdir" );
    _hist2 = tfdir.make<TH2F>( "hist2"  , "Radius of Hits vs Z of Hits", 
			       100,  0., 35., 100,  0., 45. );

    // Put a decorated message into the log file.  No endl or \n needed.
    edm::LogInfo(_messageCategory) << "Hello, world!";

    // Put an undecorated message into the log file.
    edm::LogVerbatim(_messageCategory) << "Hello, again, world!";

  }

  void Ex01InspectHits::endJob(){
    edm::LogInfo(_messageCategory) 
      << "Number of histogram entries: "
      << _hist1->GetEntries();
  }


  void Ex01InspectHits::beginRun(edm::Run const& run,
                                 edm::EventSetup const& eSetup ){

    // Access the run number.
    edm::LogInfo(_messageCategory) 
      << "Starting a new run: "
      << run.id().run();
  }

  void Ex01InspectHits::beginLuminosityBlock(edm::LuminosityBlock const& lblock,
					     edm::EventSetup const&){
    // Access the run number.
    edm::LogInfo(_messageCategory) 
      << "Starting a new lumi block: "
      << lblock.id().luminosityBlock();

  }


  void
  Ex01InspectHits::analyze(const edm::Event& evt, edm::EventSetup const&) {

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Instance name of the module that created the hits of interest;
    static const string creatorName("ex02hitmaker");

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<ToyHitCollection> handle;
    evt.getByLabel(creatorName,handle);
    
    // Some printout.
    edm::LogVerbatim log(_messageCategory);
    log << "Analyzing event #: " 
	<< evt.id().event() << ". Number of hits: "
	<< handle->size() << ".";

    // Fill histogram with number of hits per event.
    _hist4->Fill(handle->size());
    
    // Loop over all hits.
    int n(0);
    ToyHitCollection::const_iterator b = handle->begin();
    ToyHitCollection::const_iterator e = handle->end();
    for ( ; b!=e; ++b){
      
      // These are like aliases.  Used for readability.
      const ToyHit& hit = *b;
      const Hep3Vector& pos = hit._position;
      
      // Fill the histograms
      float radius = pos.perp();  
      _hist1->Fill(radius);
      _hist2->Fill(pos.z(),radius);
      _hist3->Fill(hit._pulseheight);     
      
      // Limit verbose printout.
      if ( _nAnalyzed <= _maxFullPrint ) {

	// We do need the newline here!
	log << "\nEvent: "
	    << evt.id().event() 
	    << "  | Hit #: "
	    << n++ << "  | Position: "
	    << pos
	    << " | Pulse Height: "
	    << hit._pulseheight;
      }
    }
  }
}


using mu2e::Ex01InspectHits;
DEFINE_FWK_MODULE(Ex01InspectHits);
