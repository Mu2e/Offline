//
// An EDAnalyzer Module for tuning of HoughCircles
//
// $Id: HoughTuner_module.cc,v 1.4 2011/05/18 15:47:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 15:47:40 $
//
// Original author P. Shanahan
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
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

  class HoughTuner : public art::EDAnalyzer {
  public:
    explicit HoughTuner(fhicl::ParameterSet const& pset) :
      _maxFullPrint(pset.get<int>("maxFullPrint",10)),
      _nAnalyzed(0),
      _messageCategory("ToyHitInfo"),
      _hitCreatorName(pset.get<string>("hitCreatorName"))
    { }
    virtual ~HoughTuner() { }

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run const &r);

    virtual void beginSubRun(art::SubRun const& lblock);

    // This is called for each event.
    void analyze(const art::Event& e);


  private:

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms to be filled.
    TH1F* _hBiasRadiusNoise; // mean radius of noise hits minus hough track's
    TH1F* _hBiasRadiusPhys; // mean radius of physics hits minus hough track's
    TH1F* _hBiasRadiusNoiseW; // " " wide plot
    TH1F* _hBiasRadiusPhysW; // " " wide plot

    // A category for the error logger.
    const std::string _messageCategory;

    //name of the module that created the hits to be used
    const std::string _hitCreatorName;

    void bookEventHistos(art::EventNumber_t);
    void fillEventHistos(
         mu2e::houghtransform::HoughTransform::houghCircleStruct);
  };


  void HoughTuner::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // Create some 1D histograms.
    _hBiasRadiusNoise= tfs->make<TH1F>( "hBiasRadiusNoise",
            "Mean Radius of noise hits vs. hough track",200,-20,20);
    _hBiasRadiusPhys= tfs->make<TH1F>( "hBiasRadiusPhys",
            "Mean Radius of physics hits vs. hough track",200,-20,20);
    _hBiasRadiusNoiseW= tfs->make<TH1F>( "hBiasRadiusNoiseW",
            "Mean Radius of noise hits vs. hough track (wide)",200,-1000,1000);
    _hBiasRadiusPhysW= tfs->make<TH1F>( "hBiasRadiusPhysW",
            "Mean Radius of physics hits vs. hough track (wide)",200,-1000,1000);


  }

  void HoughTuner::endJob(){
  }



  void HoughTuner::beginRun(art::Run const& run){
  }

  void HoughTuner::beginSubRun(art::SubRun const& lblock){
  }


  void HoughTuner::analyze(const art::Event& evt) {


    static int ncalls(0);
    ++ncalls;

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Ask the event to give us a "handle" to the requested hits.
    //    art::Handle<StepPointMCCollection> hits;
    //evt.getByLabel(creatorName,hits);

    static const string collectionName("tracker");
    art::Handle<StepPointMCCollection> hitsHandle;
    evt.getByLabel(_hitCreatorName,collectionName,hitsHandle);
    StepPointMCCollection const* hits = hitsHandle.product();

    art::Handle<HoughCircleCollection> hcHandle;
    evt.getByType(hcHandle);

   // primary hough circle: i.e., the first one in the collection
    const HoughCircle* hcp=0;
    if (hcHandle->size()) hcp=&hcHandle->at(0);

    // Fill histogram with number of hits per event.
    //_hHitMultiplicity->Fill(hits->size());

    // and number of Hough Tracks
   // _hHoughMultiplicity->Fill(hcHandle->size());

    // and make a pretty Plot
     unsigned evtno=evt.id().event();

     RootNameTitleHelper displayCanv("canvEvt","Display, Event",evtno,5);
     art::ServiceHandle<art::TFileService> tfs;

     TCanvas *canvEvt = tfs->make<TCanvas>(displayCanv.name(),
                   displayCanv.title());
     TH2F hpad("event","Event Display",100,-1000.,1000.,100,-1000,1000);
     canvEvt->cd();
     hpad.Draw();

// draw circles
     std::vector<TEllipse*> circles; // to be able to delete them when done
std::cout<<"nCircles="<<hcHandle->size()<<std::endl;
     for (size_t ihc=0; ihc<hcHandle->size(); ihc++) {
         const HoughCircle& hc=hcHandle->at(ihc);
         TEllipse *circle=new TEllipse(hc.Center().x(),hc.Center().y(),hc.Radius());
         circles.push_back(circle);
         circle->Draw();
     }


// draw hits
     enum {kPhys,kNoise};
     TPolyMarker* tp[2]={0,0};
     tp[kPhys]=new TPolyMarker(); tp[kPhys]->SetMarkerStyle(20);
            tp[kPhys]->SetMarkerColor(kRed);
     tp[kNoise]=new TPolyMarker(); tp[kNoise]->SetMarkerStyle(24);
            tp[kNoise]->SetMarkerColor(kCyan);
     int npm[2]={0,0};

     // also want to calculate mean hit radius from Hough Center
     double rmean[2]={0,0};

// loop over hits
     for (size_t ih=0; ih<hits->size(); ih++)
     {

         const StepPointMC& hit=hits->at(ih);


          // add hit to noise or physics TPolyMarker and stats as appropriate
         unsigned int kType=kPhys;
         if ( 2==hit.trackId().asInt() && TMath::Abs(hit.totalEDep()-0)<1e-10)
                                                              kType=kNoise;
         tp[kType]->SetPoint(npm[kType]++,hit.position()[0],hit.position()[1]);

         // calculate mean hit radius from center of 1st circle
         if (hcp) {
           double rhit=sqrt(TMath::Power(hit.position()[0]-hcp->Center().x(),2)
                           +TMath::Power(hit.position()[1]-hcp->Center().y(),2));
           rmean[kType]+=rhit;
         }

     }

     // normalize the means
     if(npm[kNoise]) rmean[kNoise]/=npm[kNoise];
     if(npm[kPhys]) rmean[kPhys]/=npm[kPhys];

     tp[kPhys]->Draw();
     tp[kNoise]->Draw();
     canvEvt->Write();

     delete tp[kPhys];
     delete tp[kNoise];
     for (size_t ic=0; ic<circles.size(); ic++) delete circles[ic];


     if (hcp) {
       _hBiasRadiusNoise->Fill(rmean[kNoise]-hcp->Radius());
       _hBiasRadiusPhys->Fill(rmean[kPhys]-hcp->Radius());
       _hBiasRadiusNoiseW->Fill(rmean[kNoise]-hcp->Radius());
       _hBiasRadiusPhysW->Fill(rmean[kPhys]-hcp->Radius());
      if (TMath::Abs(rmean[kPhys]-hcp->Radius())>100)
       std::cout<<"Event "<< evtno <<"diff="<<rmean[kPhys]-hcp->Radius()<<std::endl;
     }



  } // end of ::analyze.


  void HoughTuner::bookEventHistos(art::EventNumber_t evtno)
  {
  } //bookEventHistos

 void HoughTuner::fillEventHistos(
         mu2e::houghtransform::HoughTransform::houghCircleStruct houghCircle)
 {

 } // fillEventHistos()


} // namespace HoughTuner


using mu2e::HoughTuner;
DEFINE_ART_MODULE(HoughTuner);
