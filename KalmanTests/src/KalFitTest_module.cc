//
// Module to perform BaBar Kalman fit
//
// $Id: KalFitTest_module.cc,v 1.7 2011/06/17 21:56:57 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/17 21:56:57 $
//

// framework
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/DetStrawHitElem.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"

// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>


using namespace std; 

namespace mu2e 
{
  class KalFitTest : public art::EDProducer
  {
  public:
    explicit KalFitTest(fhicl::ParameterSet const&);
    virtual ~KalFitTest();
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    void endJob();
    void produce(art::Event& e);

  private:
    // configuration parameters
    int _diag,_debug;
    int _printfreq;
    // event object labels
    std::string _strawhitslabel;
    std::string _mcstrawhitslabel;
    std::string _mcptrlabel;
    std::string _mcstepslabel;
    // cache of event objects
    const StrawHitCollection* _strawhits;
    // mc data
    MCEvtData _mcdata;
    // fitter
    KalFit _kfit;
    // helper functions
    bool findData(art::Event& e);
    bool findMC(art::Event& e);
    // MC tools
    KalFitMC _kfitmc;
// 
  };
  
  KalFitTest::KalFitTest(fhicl::ParameterSet const& pset) : 
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",10)),
    _strawhitslabel(pset.get<std::string>("strawHitsLabel","makeSH")),
    _mcstrawhitslabel(pset.get<std::string>("MCStrawHitsLabel","makeSH")),
    _mcptrlabel(pset.get<std::string>("MCPtrLabel","makeSH")),
    _mcstepslabel(pset.get<std::string>("MCStepsLabel","g4run")),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC"))
  {
  }

  KalFitTest::~KalFitTest(){}

  void KalFitTest::beginJob(){}

  void KalFitTest::beginRun(art::Run& ){}

  void KalFitTest::produce(art::Event& event ) 
  {
// event printout
    int iev=event.id().event();
    if((iev%_printfreq)==0)cout<<"KalFitTest: event="<<iev<<endl;
// find the data
    if(!findData(event)){
      cout << "No straw hits found " << endl;
      return;
    }
// find mc truth
    if(!findMC(event)){
      cout<<"MC information missing "<< endl;
      return;
    }
    TrkDef mytrk(_strawhits);
// must initialize t0, momentum, initial trajectory.  These should come from patrec
// that doesn't yet exist. For now, take from the MC truth.  There must be a better way to define the primary particle, FIXME!!!!!!
    cet::map_vector_key trkid(1);
    if(_kfitmc.trkFromMC(_mcdata,trkid,mytrk)){
// use this to create a track
      TrkKalFit myfit;
      _kfit.makeTrack(mytrk,myfit);
// test if fit succeeded
      if(myfit._fit.success()){
//  diagnostics
        if(_diag > 0){
          _kfitmc.trkDiag(_mcdata,mytrk,myfit);
          if(_diag > 1){
            for(std::vector<TrkStrawHit*>::iterator ihit=myfit._hits.begin();ihit!=myfit._hits.end();ihit++){
              TrkStrawHit* trkhit = *ihit;
              _kfitmc.hitDiag(_mcdata,trkhit);
            }
          }
        }
      }
// cleanup; the track should be put in the event
      myfit.deleteTrack();
    }
  }
  
  void KalFitTest::endJob()
  {
// does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }
  
// find the input data objects needed to make tracks
  bool KalFitTest::findData(art::Event& evt){
    _strawhits = 0;
    art::Handle<mu2e::StrawHitCollection> strawhitsH;
    if(evt.getByLabel(_strawhitslabel,strawhitsH))
      _strawhits = strawhitsH.product();
    return _strawhits != 0;
  }
  
// find the MC truth objects in the event and set the local cache
  bool KalFitTest::findMC(art::Event& evt) {
    _mcdata.clear();
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    if(evt.getByLabel(_mcstrawhitslabel,truthHandle))
      _mcdata._mcstrawhits = truthHandle.product();
  // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> mchitptrHandle;
    if(evt.getByLabel(_mcptrlabel,"StrawHitMCPtr",mchitptrHandle))
      _mcdata._mchitptr = mchitptrHandle.product();
  // Get the persistent data about the StepPointMCs, from the tracker and the virtual detectors
    art::Handle<StepPointMCCollection> mctrackerstepsHandle;
    if(evt.getByLabel(_mcstepslabel,"tracker",mctrackerstepsHandle))
      _mcdata._mcsteps = mctrackerstepsHandle.product();
    art::Handle<StepPointMCCollection> mcVDstepsHandle;
    if(evt.getByLabel(_mcstepslabel,"virtualdetector",mcVDstepsHandle))
      _mcdata._mcvdsteps = mcVDstepsHandle.product();
    return _mcdata.good();
  }
}

using mu2e::KalFitTest;
DEFINE_ART_MODULE(KalFitTest);
