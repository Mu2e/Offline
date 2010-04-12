//
// An EDAnalyzer Module for tuning of HoughCircles
//
// $Id: HoughTuner_plugin.cc,v 1.2 2010/04/12 18:26:12 shanahan Exp $
// $Author: shanahan $ 
// $Date: 2010/04/12 18:26:12 $
//
// Original author P. Shanahan
//

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
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TPolyMarker.h"
#include "TMath.h"
#include "TF1.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"


// Mu2e includes.
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/HoughCircleCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "HoughTransform/inc/HoughTransform.hh"
#include "HitCluster/inc/HitCluster.hh"
#include "GeneralUtilities/inc/RootNameTitleHelper.hh"
#include "GeneralUtilities/inc/pow.hh"

//CLHEP includes
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Randomize.h"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::RandPoisson;
using namespace mu2e;
using namespace mu2e::houghtransform;
namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 

  class HoughTuner : public edm::EDAnalyzer {
  public:
    explicit HoughTuner(edm::ParameterSet const& pset) : 
      _hitCreatorName(pset.getParameter<string>("hitCreatorName")),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",10)),
      _nAnalyzed(0),
      _messageCategory("ToyHitInfo")
    { }
    virtual ~HoughTuner() { }

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
    TH1F* _hRadius;

    // A category for the error logger.
    const std::string _messageCategory;

    //name of the module that created the hits to be used
    const std::string _hitCreatorName;

    void HoughTuner::bookEventHistos(edm::EventNumber_t);
    void HoughTuner::fillEventHistos(
         mu2e::houghtransform::HoughTransform::houghCircleStruct);
  };


  void HoughTuner::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius", "Radius of Hits;(mm)",          100,  0., 1000. );

  }

  void HoughTuner::endJob(){
  }



  void HoughTuner::beginRun(edm::Run const& run,
                                 edm::EventSetup const& eSetup ){
  }

  void HoughTuner::beginLuminosityBlock(edm::LuminosityBlock const& lblock,
					     edm::EventSetup const&){
  }


  void HoughTuner::analyze(const edm::Event& evt, edm::EventSetup const&) {


    static int ncalls(0);
    ++ncalls;

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Ask the event to give us a "handle" to the requested hits.
    //    edm::Handle<StepPointMCCollection> hits;
    //evt.getByLabel(creatorName,hits);
    
    edm::Handle<StepPointMCCollection> hitsHandle;
    evt.getByLabel(_hitCreatorName,hitsHandle);
    StepPointMCCollection const* hits = hitsHandle.product();

    edm::Handle<HoughCircleCollection> hcHandle;
    evt.getByType(hcHandle);

    // Fill histogram with number of hits per event.
    //_hHitMultiplicity->Fill(hits->size());
 
    // and number of Hough Tracks
   // _hHoughMultiplicity->Fill(hcHandle->size());

    // and make a pretty Plot
     uint32_t evtno=evt.id().event();

     RootNameTitleHelper displayCanv("_canvEvt","Display, Event",evtno,5);
     edm::Service<edm::TFileService> tfs;

     TCanvas *canvEvt = tfs->make<TCanvas>(displayCanv.name(),
                   displayCanv.title());
     TH2F hpad("event","Event Display",100,-1000.,1000.,100,-1000,1000);
     canvEvt->cd();
     hpad.Draw();

// draw circles
     std::vector<TEllipse*> circles; // to be able to delete them when done
     for (int ihc=0; ihc<hcHandle->size(); ihc++) {
         const HoughCircle& hc=hcHandle->at(ihc);
         TEllipse *circle=new TEllipse(hc.Center().x(),hc.Center().y(),hc.Radius());
         circles.push_back(circle);
         circle->Draw();
        std::cout<<"Draw a circle: "<<hc.Center().x()<<" "<<hc.Center().y()
                <<" "<<hc.Radius()<<std::endl;
     }


// draw hits
     TPolyMarker* tpPhys=new TPolyMarker(); tpPhys->SetMarkerStyle(20);
            tpPhys->SetMarkerColor(kRed);
     TPolyMarker* tpNois=new TPolyMarker(); tpNois->SetMarkerStyle(24);
            tpNois->SetMarkerColor(kCyan);
     int nPhys=0;
     int nNois=0;
     
     for (int ih=0; ih<hits->size(); ih++)
     {
          
         const StepPointMC& hit=hits->at(ih);

          // add hit to noise or physics TPolyMarker as appropriate
         if ( 2==hit.trackId() && TMath::Abs(hit.totalEDep()-0)<1e-10) {
             tpNois->SetPoint(nNois++,hit.position()[0],hit.position()[1]);
         } else {
             tpPhys->SetPoint(nPhys++,hit.position()[0],hit.position()[1]);
         }

     }

     tpPhys->Draw();
     tpNois->Draw();
     canvEvt->Write();

     delete tpPhys;
     delete tpNois;
     for (int ic=0; ic<circles.size(); ic++) delete circles[ic];
     

 
  } // end of ::analyze.


  void HoughTuner::bookEventHistos(edm::EventNumber_t evtno)
  {	
  } //bookEventHistos
  
 void HoughTuner::fillEventHistos(
         mu2e::houghtransform::HoughTransform::houghCircleStruct houghCircle) 
 {
  
 } // fillEventHistos()


} // namespace HoughTuner


using mu2e::HoughTuner;
DEFINE_FWK_MODULE(HoughTuner);
