//
// Module which starts the event display, and transmits the data of each event to the event display.
//
// $Id: kalmanFitv2_module.cc,v 1.4 2012/12/04 00:51:27 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2012/12/04 00:51:27 $
//

#include <iostream>
#include <string>
#include <memory>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "ITrackerGeom/inc/ITracker.hh"
//#include "GeometryService/inc/WorldG4.hh"
//#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"

#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"

#include <iostream>
#include <fstream>

#include <iostream>
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/KalFitResult.hh"
//#include "KalmanTests/inc/TrkDef.hh"
//#include "KalmanTests/inc/TrkStrawHit.hh"
//#include "KalmanTests/inc/KalFit.hh"
//#include "KalmanTests/inc/KalFitMC.hh"
//#include "TrkPatRec/inc/TrkHitFilter.hh"
//#include "TrkPatRec/inc/TrkHelixFit.hh"
//#include "TrkBase/TrkPoca.hh"
//#include "TrkPatRec/src/TrkPatRec_module.cc"
#include "TrkBase/HelixParams.hh"
#include "TrkPatRec/inc/TrkPatRec.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkPoca.hh"
// Data Output
#include "KalmanTests/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"

#include "DchGeom/DchDetector.hh"

#include "KalmanTestsI/inc/DchGDchP.hh"
#include "KalmanTestsI/inc/KalFitI.hh"

#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <fstream>

#include "CLHEP/Random/RandGaussQ.h"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "KalmanTestsI/inc/kalFitOutUtils.hh"

using namespace std; 

namespace mu2e 
{
  class kalmanFitv2 : public art::EDProducer
  {
    public:
    explicit kalmanFitv2(fhicl::ParameterSet const&);
    virtual ~kalmanFitv2() { }
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    void endJob();
    void produce(art::Event& e);

    private:
    CLHEP::RandGaussQ _gaussian;


    // Label of the module that created the data products.
    std::string _g4ModuleLabel;
    // Label of the generator.
    std::string _generatorModuleLabel;
    // Instance names of data products
    std::string _targetStepPoints;


    string _hitMakerModuleLabel;

    
    double _cutdp;
    double _aveSpcRes;
    bool _usemchitdist;
    
    int _debugLvl;
    bool _writeDataOut;

    DchDetector* dchdet;
    std::string _materialdb;
    bool _doEndCapExtrapol;
    TrkParticle _tpart; // particle type being searched for
    TrkFitDirection _fdir;  // fit direction in search

    KalFitI _kfit;

    std::string _iname; // data instance name
    //
    PayloadSaver _payloadSaver;
    // End: run time parameters

//    bool FillMCInfo(art::Event& event);

    std::vector<hitIndex> strawhits;
    HelixTraj seed;
//    double time0,time0_max;
//    double s0;
//    std::vector<double> hitflt;

//    CLHEP::Hep3Vector gen_pos;
//    CLHEP::Hep3Vector gen_mom;

    kalFitOutUtils _histoOut;

    bool _firstEv;
    Int_t _eventid;

  };

  kalmanFitv2::kalmanFitv2(fhicl::ParameterSet const& pset) :
    _gaussian( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _targetStepPoints(pset.get<string>("targetStepPoints","stoppingtarget")),
    _hitMakerModuleLabel(pset.get<std::string>("hitMakerModuleLabel", "makeDcH")),
    _cutdp(pset.get<double>("cutdp",1000.)),
    _aveSpcRes(pset.get<double>("aveSpcRes",0.2)),
    _usemchitdist(pset.get<bool>("mchitdist",false)),
    _debugLvl(pset.get<int>("debuglevel",1)),
    _writeDataOut(pset.get<bool>("writeDataOut",false)),
    _materialdb(pset.get<std::string>("materialdb", "")),
    _doEndCapExtrapol(pset.get<bool>("doEndCapExtrapol",false)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit")),
    _payloadSaver(pset),
    seed(TrkParams(HelixTraj::NHLXPRM)),
    _histoOut (_g4ModuleLabel,_generatorModuleLabel,_targetStepPoints,_hitMakerModuleLabel,
               _debugLvl,pset.get<int>("oneturn",0),pset.get<string>("vdStepPoints","virtualdetector"),
               GenId(pset.get<int>("genidselect",GenId::conversionGun))),
    _firstEv(true)
//  ,
//    nbadfit(0),
//    startpar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1)),
//    recopar(CLHEP::HepVector(5,0),CLHEP::HepSymMatrix(5,1))
  {
    produces<CLHEP::HepVector>("test");

    _iname = _fdir.name() + _tpart.name();
    if (_writeDataOut) {
            produces<KalRepCollection>(_iname);
            produces<KalRepPayloadCollection>();
    }
// set # bins for time spectrum plot
  }

  void kalmanFitv2::beginJob()
  {
         _histoOut.bookHitos();
         _eventid = 0;
  }

  void kalmanFitv2::beginRun(art::Run&){
    std::cout<<"init dchdetector"<<std::endl;
    DchGDchP gdch(_materialdb,_kfit.useDetailedWrSuppDescr(),_kfit.useSingleCellDescr());
    dchdet=new DchDetector(gdch,true);
  }

  void kalmanFitv2::produce(art::Event& event) 
  {
    if (_firstEv) {_histoOut._bfield=&_kfit.bField(); _firstEv=false;}
    ++_eventid;

    if(_debugLvl>9) 
      ErrMsg::ErrLoggingLevel = debugging;


    int iev=event.id().event();
    if (_debugLvl>9) { cout<<"event="<<iev<<endl; }
    else if((iev%500)==0) { cout<<"event="<<iev<<endl; }

    art::Handle<mu2e::StrawHitMCTruthCollection> cellmchits;
    art::Handle<mu2e::StrawHitCollection> cellhits;
    if(!(event.getByLabel(_hitMakerModuleLabel,cellhits)&&event.getByLabel(_hitMakerModuleLabel,cellmchits))){
            if (_debugLvl>0) {
                    cout<<"no dc hits "
                                    <<event.getByLabel(_hitMakerModuleLabel,cellhits)<<" "
                                    <<event.getByLabel(_hitMakerModuleLabel,cellmchits)<<endl;
            }
      return;
    }

    if(!_histoOut.FillMCInfo(event,strawhits,seed)) return;
    
    // create output
    auto_ptr<KalRepCollection> tracks(new KalRepCollection );

    TrkDef trkdef(cellhits.product(),strawhits,seed,_tpart,_fdir);//,TrkParticle(),TrkFitDirection());
    trkdef.setT0(TrkT0(_histoOut.recoinfo.t0,0.0));

    _kfit._flt0=_histoOut.s0;
    _kfit._hitflt=_histoOut.hitflt;
    KalFitResult myfit(trkdef);
    _kfit.makeTrack(myfit);

    if(myfit._fit.success()){
     TrkHotList* hots = myfit._krep->hotList();

     if (_debugLvl>0) { std::cout<<"active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }
     //std::vector<hitIndex> strawhitsinactive;
      
     _kfit.reActivateHitsbyTurn(myfit);
      
     if (_debugLvl>9) { std::cout<<"ReTurn:active hits "<<hots->nActive()<<" nhit "<<hots->nHit()<<std::endl; }
    }


    if(_doEndCapExtrapol && myfit._fit.success()) {
            _kfit.makeExtrapolOutTrack(myfit);
    }

    if(myfit._fit.success()){
      auto_ptr<HepVector> recomom(new HepVector(6));

      //KalRep *kalrep = myfit._krep;
      double fltlen=myfit._krep->firstHit()->globalLength();
      HepVector hvpos(3),hvmom(3);
      HepSymMatrix ehpos(3),ehmom(3);
      HepMatrix ehposmom(3,3);
      myfit._krep->getAllWeights(fltlen,hvpos,hvmom,ehpos,ehmom,ehposmom);
      for(int i=1;i<=3;i++) {
	(*recomom)(i)=hvpos(i);
	(*recomom)(i+3)=hvmom(i);
      }
    
      if (_writeDataOut) {
        event.put(recomom,"test");
      }
    }

    _histoOut.FillHistos(myfit,seed);

    if(myfit._fit.success()){
            // save successful kalman fits in the event
            tracks->push_back( myfit.stealTrack() );
    } else {
            myfit.deleteTrack();
    }

    if (_writeDataOut) {
            // put the tracks into the event
            art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
            _payloadSaver.put(*tracks, tracksID, event);
            event.put(tracks,_iname);
    }

  }

  void kalmanFitv2::endJob()
  {
          _histoOut.finalizeHistos();
  }

}  // end namespace mu2e

using mu2e::kalmanFitv2;
DEFINE_ART_MODULE(kalmanFitv2);
