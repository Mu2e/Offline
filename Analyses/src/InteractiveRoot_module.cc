//
// A plugin to show how to use interactive ROOT with the framework.
//
// $Id: InteractiveRoot_module.cc,v 1.9 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Rob Kutschke
//
// Notes
// 1) The key to making this work is the management of the TApplication instance.
//    In order for ROOT interactive graphics to work, exactly one instance
//    of TApplication must be created somewhere in the job.  If several
//    modules wish to do ROOT graphics then exactly one of them needs to
//    create the TApplication; which ever module does this must also
//    delete the TApplication at end of job.  In this module, the test
//    to see if the TApplication was created relies on the value of gApplication.
//    The unique_ptr will look after deletion at destructor-time.
//
// 2) ROOT requires that TCanvas objects have unique names, even if we create
//    them in different modules or write them to different TDirectory's. If we
//    don't respect this, ROOT will crash.  I presume this is because TCanvas's
//    are all owned by the TApplication, which no knowledge of modules or TDirectories.
//    In the framework the module label is guaranteed to be unique within a job.
//    So this example embeds the module label into the canvas name.  Making the titles
//    unique if for convenience.
//
// 3) There is a deficiency in TFileService. It correctly does management of ROOT directories
//    for beginJob, analyze and some other methods. It does not do this correctly for
//    endJob.  So we need to do the directory management ourselves, hence the existence
//    of the _directory data member.  This must be initialized after TFileService has
//    created the directory for this module, which occurs on the first call to tfs->make<T>(...).
//    This should be corrected in a new version of TFileService.
//
// 4) For background see http://mu2e.fnal.gov/atwork/computing/ROOTFAQ.shtml .
//    There is something broken in the interactions among ROOT, TFileService, the EDM and the
//    framework.  I expected to be able to code the following right after instantiating the
//    canvas:
//       gDirectory->Append(_canvas)
//    Then we would not need anything in the endJob method.  I am not sure why but
//    this causes a crash at the end of the job; the current guess is that there is a
//    double delete of the histogram in the canvas, once when the histogram is deleted and
//    once when the canvas is deleted.  The hack solution is to do a Write() in endJob.
//

// C++ includes.
#include <iostream>
#include <string>
#include <memory>
#include <unistd.h>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class Straw;

  class InteractiveRoot : public art::EDAnalyzer {
  public:

    explicit InteractiveRoot(fhicl::ParameterSet const& pset);
    virtual ~InteractiveRoot() { }

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(const art::Event& e);

  private:

    // Start: run time parameters

    // The module label of this module.
    std::string _moduleLabel;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Number of events to accumulate between prompts.
    int _nAccumulate;

    // End: run time parameters

    // Pointers to histograms, ntuples, TGraphs.
    TH1F*         _hMultiplicity;
    TCanvas*      _canvas;

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  InteractiveRoot::InteractiveRoot(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _moduleLabel(pset.get<string>("module_label")),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _nAccumulate(pset.get<int>("nAccumulate",20)),

    // ROOT objects that are the main focus of this example.
    _hMultiplicity(0),
    _canvas(0),

    // Some ugly but necessary ROOT related bookkeeping.
    _application(nullptr),
    _directory(0){

  }

  void InteractiveRoot::beginJob( ){

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    // Create a histogram.
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Hits per Event", 100,  0.,  100. );

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
    }

    // Create a canvas with a unique name.  See note 2.
    TString name  = "canvas_"     + _moduleLabel;
    TString title = "Canvas for " + _moduleLabel;
    int window_size(700);
    _canvas = tfs->make<TCanvas>(name,title,window_size,window_size);

    // Draw the still empty histogram. It will be updated later.
    _hMultiplicity->Draw();

    // See note 3.
    _directory = gDirectory;

  }

  void InteractiveRoot::analyze(const art::Event& event) {

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hitsHandle;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hitsHandle);
    StepPointMCCollection const& hits = *hitsHandle;

    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits.size());

    // Periodically update the displayed histogram.
    if ( event.id().event()%_nAccumulate==0 ){
      _canvas->Modified();
      _canvas->Update();

      cerr << "Sleeping 5 seconds for you to have a look at the canvas: "
           << _moduleLabel << endl;
      sleep(5);

    }

  } // end analyze

  void InteractiveRoot::endJob(){

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }

}  // end namespace mu2e

using mu2e::InteractiveRoot;
DEFINE_ART_MODULE(InteractiveRoot);
