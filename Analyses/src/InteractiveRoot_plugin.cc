//
// A plugin to test using root interactively.
//
// $Id: InteractiveRoot_plugin.cc,v 1.1 2010/09/15 23:14:13 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/09/15 23:14:13 $
//
// Original author Rob Kutschke
//


// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class Straw;

  class InteractiveRoot : public edm::EDAnalyzer {
  public:
    
    explicit InteractiveRoot(edm::ParameterSet const& pset);
    virtual ~InteractiveRoot() { }

    virtual void beginJob(edm::EventSetup const&);
    void endJob();
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Start: run time parameters

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Pointers to histograms, ntuples, TGraphs.
    TH1F*         _hMultiplicity;
    TCanvas*      _canvas;
    TApplication* _application;
    
  };

  InteractiveRoot::InteractiveRoot(edm::ParameterSet const& pset) : 

    // Run time parameters
    _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),

    // ROOT objects
    _hMultiplicity(0),
    _application(0),
    _canvas(0){
  }
  
  void InteractiveRoot::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    // Create a histogram.
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Hits per Event", 100,  0.,  100. );

    // Create a root interactive environment.
    // This may be done only once per job.  If multiple modules need to make TCanvases, then
    // we need to find a way to coordinate this step; move it to the TFileService?
    int    tmp_argc(0);
    char** tmp_argv(0);
    _application = new TApplication( "noapplication", &tmp_argc, tmp_argv );

    // Create a canvas
    _canvas = new TCanvas("canvas","My Canvas",200,10,700,700);

    // Draw the still empty histogram.
    _hMultiplicity->Draw();

  }

  void InteractiveRoot::analyze(const edm::Event& event, edm::EventSetup const&) {

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hitsHandle;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hitsHandle);
    StepPointMCCollection const& hits = *hitsHandle;

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits.size());

    // Periodically update the displayed histogram.
    if ( event.id().event()%10==0 ){
      _canvas->Modified();
      _canvas->Update();

      cerr << "Double click in the Canvas to continue:" ;
      _canvas->WaitPrimitive();
      cerr << endl;

    }

  } // end analyze

  void InteractiveRoot::endJob(){

  }

}  // end namespace mu2e

using mu2e::InteractiveRoot;
DEFINE_FWK_MODULE(InteractiveRoot);
