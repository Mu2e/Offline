//
// Module to perform BaBar Kalman fit
//
// $Id: KalFitTest_module.cc,v 1.15 2012/07/23 22:30:57 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/23 22:30:57 $
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
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
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
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
   // cache of event objects
    const StrawHitCollection* _strawhits;
    // fitter
    KalFit _kfit;
    // helper functions
    bool findData(art::Event& e);
    // MC tools
    KalFitMC _kfitmc;
// 
  };
  
  KalFitTest::KalFitTest(fhicl::ParameterSet const& pset) : 
    _diag(pset.get<int>("diagLevel",0)),
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",10)),
    _strawhitslabel(pset.get<std::string>("strawHitsLabel","makeSH")),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC"))
  {
    produces<KalRepCollection>();
  }

  KalFitTest::~KalFitTest(){}

  void KalFitTest::beginJob(){
     if(_diag > 0)
      TTree* trkdiag = _kfitmc.createTrkDiag();
  }

  void KalFitTest::beginRun(art::Run& ){
  }

  void KalFitTest::produce(art::Event& event ) 
  {
    auto_ptr<KalRepCollection> tracks(new KalRepCollection );
// event printout
    int iev=event.id().event();
    if((iev%_printfreq)==0)cout<<"KalFitTest: event="<<iev<<endl;
// find the data
    if(!findData(event)){
      cout << "No straw hits found " << endl;
      return;
    }
// find mc truth
    if(!_kfitmc.findMCData(event)){
      cout<<"MC information missing "<< endl;
      return;
    }
    TrkDef mytrk(_strawhits);
    // must initialize t0, momentum, initial trajectory.  These should come from patrec
    // that doesn't yet exist. For now, take from the MC truth.  There must be a better way to define the primary particle, FIXME!!!!!!
    cet::map_vector_key trkid(1);
    if(_kfitmc.trkFromMC(trkid,mytrk)){
      // use this to create a track
      TrkKalFit myfit;
      _kfit.makeTrack(mytrk,myfit);
      //  diagnostics
      if(_diag > 0){
	_kfitmc.kalDiag(myfit._krep);
	if(_diag > 1){
	  for(std::vector<TrkStrawHit*>::iterator ihit=myfit._hits.begin();ihit!=myfit._hits.end();ihit++){
	    TrkStrawHit* trkhit = *ihit;
	    _kfitmc.hitDiag(trkhit);
	  }
	}
      }
      // If fit is successful, pass ownership of the track to the event.
      if(myfit._krep != 0)tracks->push_back( myfit.stealTrack() );
    }
    event.put(tracks);
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
  
}

using mu2e::KalFitTest;
DEFINE_ART_MODULE(KalFitTest);
