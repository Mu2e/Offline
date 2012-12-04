//
// Patter recognition for the ITracker (based on ILC 4th PR)
//
// $Id: ITTrackMiniIReco_module.cc,v 1.1 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
//

// C++ includes.
#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "ITrackerGeom/inc/Cell.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "TrkPatRec/inc/TrkPatRec.hh"

#include "TrkPatRecI/inc/fillWires.hh"
#include "TrkPatRecI/miniilc/src/IlcDCHtracker.h"
#include "TrkPatRecI/miniilc/src/IlcDCHseed.h"
#include "TrkPatRecI/miniilc/src/IlcDCHParam.h"
#include "TrkPatRecI/miniilc/src/IlcDCHReconstructor.h"
#include "TrkPatRecI/miniilc/src/IlcLog.h"
#include "TrkPatRecI/miniilc/src/IlcPID.h"
#include "TrkPatRecI/miniilc/src/readData.h"

// ROOT includes
#include "TRandom.h"

#define invSqrt12 0.2886751345948129

namespace mu2e {

  class ITTrackMiniIReco : public art::EDProducer{
  public:
    
    explicit ITTrackMiniIReco(fhicl::ParameterSet const& pset);
    virtual ~ITTrackMiniIReco() {};
    
    virtual void beginJob();
    void beginRun();
    int firstEvent;
    void endJob();
    void produce(art::Event & e);

  private:

    // Start: run time parameters
     IlcDCHtracker *itracker;

     int _debugLevel;
    
    // The module labels
    std::string _moduleLabel;
    std::string _hitmakerModuleLabel;
    std::string _timeRejecterModuleLabel;
    bool _useZCoordinate;

    bool _doMinPrints;
  };

  ITTrackMiniIReco::ITTrackMiniIReco(fhicl::ParameterSet const& pset) :
    // Run time parameters
    _debugLevel(pset.get<int>("debugLevel",0)),
    _moduleLabel(pset.get<std::string>("module_label")),
    _hitmakerModuleLabel(pset.get<std::string>("hitMakerModuleLabel", "makeDcH")),
    _timeRejecterModuleLabel(pset.get<std::string>("tRejecterModuleLabel","trckRjctByTime"))  ,
    _useZCoordinate(pset.get<bool>("useZCoordinate",false))
  {
    // Tell the framework what we make.
    produces<TrackSeedCollection>();
    firstEvent=0;
    _doMinPrints = (_debugLevel>0);
  }
  
  void ITTrackMiniIReco::beginJob(){
    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;
  }
 
  void ITTrackMiniIReco::beginRun(){

    // Get access to the TFile service and save current directory for later use.

    IlcTracker::SetBz(-10.);
    new IlcDCHReconstructor;//initialized reco params
    IlcDCHParam *par=new IlcDCHParam;
    //IlcPID *pid=new IlcPID;//initialize particles
    //pid->GetName();
    par->SetExperiment(1);
    par->SetEndCapType(0);
    IlcDCHwireposition *wpos= fillWires();
    itracker=new IlcDCHtracker(par,wpos);
    itracker->SetDebug(_debugLevel);
    if(_useZCoordinate)itracker->SetUseZCoordinate();
    //    IlcLog::SetGlobalDebugLevel(20);

  }

  void ITTrackMiniIReco::produce(art::Event & event ) {
    if(!firstEvent){firstEvent=1;beginRun();}

    
    auto_ptr<TrackSeedCollection> outseeds(new TrackSeedCollection);


    //const Tracker& tracker = getTrackerOrThrow();
    //const ITracker &itr = static_cast<const ITracker&>( tracker );
    //CellGeometryHandle *itwp = itr.getCellGeometryHandle();
    

    art::Handle<mu2e::PtrStepPointMCVectorCollection> cellptr;
    event.getByLabel(_hitmakerModuleLabel,"StrawHitMCPtr",cellptr);
    PtrStepPointMCVectorCollection const* hits_mcptr = cellptr.product();
    if(hits_mcptr->size()==0) {
      if (_doMinPrints) { std::cout<<"Empty ev = "<<event.id().event()<<std::endl; }
      return;
    }
    
    art::Handle<mu2e::StrawHitMCTruthCollection> cellmchits;
    art::Handle<mu2e::StrawHitCollection> cellhits;
    if(!(event.getByLabel(_hitmakerModuleLabel,cellhits)&&event.getByLabel(_hitmakerModuleLabel,cellmchits))){
            if (_doMinPrints) {
                    std::cout<<"no dc hits "
                                    <<event.getByLabel(_hitmakerModuleLabel,cellhits)<<" "
                                    <<event.getByLabel(_hitmakerModuleLabel,cellmchits)<<std::endl;
            }
      return;
    }
    StrawHitCollection const* shits = cellhits.product();
    StrawHitMCTruthCollection const* shitmcs = cellmchits.product();
    
    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();

    int firstmcstep=-1;
    double time0_min=1e10;
    double time0_max=-1e10;
    for(unsigned int i=0;i<hits_mcptr->size();i++){
      const PtrStepPointMCVector&  mcptr=hits_mcptr->at(i);
      double tcell=mcptr.at(0)->time();
      if(mcptr.size()&&mcptr.at(0)->trackId().asInt()==1&&
	 tcell<time0_min){time0_min=tcell;firstmcstep=i;}
      if(mcptr.size()&&mcptr.at(0)->trackId().asInt()==1&&
	 tcell>time0_max){time0_max=tcell;}
    }
    if(firstmcstep<0){
            if (_doMinPrints) {
                    std::cout<<"ev = "<<event.id().event()
                                          <<" Not electron in all cells, ncell="<<hits_mcptr->size();
                    if(hits_mcptr->size()&&hits_mcptr->at(0).size()) {
                            std::cout<<" id= "<< hits_mcptr->at(0).at(0)->trackId().asInt();
                    }
                    std::cout<<std::endl;
            }
      return;
    }

    const PtrStepPointMCVector& mcptr=hits_mcptr->at(firstmcstep);
    double  time0=mcptr.at(0)->time();



    clock_t startClock = clock();

    size_t nStrawPerEvent = shits->size();
    size_t nTimeClusPerEvent = tclusts->size();
    if (_doMinPrints) {
            std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
            std::cout<<"event "<<event.id().event()<<" tot N hit "<<nStrawPerEvent<<" N peak found "<<nTimeClusPerEvent<<std::endl;
            std::cout<<"----------------------------------------------------------------------------------"<<std::endl;
    }
    
    std::vector<TrkTimePeak> tpeaks;
    
    for (size_t ipeak=0; ipeak<nTimeClusPerEvent; ipeak++) {
      
      TrackerHitTimeCluster const&  tclust(tclusts->at(ipeak));

      double tpeakerr = (tclust._maxHitTime - tclust._minHitTime)*invSqrt12;
      if (_doMinPrints) { std::cout<<ipeak<<" nhits in peak "<<tclust._selectedTrackerHits.size()<<" from "<<shits->size()
	       <<" t= "<<tclust._meanTime<<" +- "<<tpeakerr<<" min "<<tclust._minHitTime<<" max "<<tclust._maxHitTime
	       <<" t0mc "<<time0<<" t0maxmc "<<time0_max<<std::endl;
      }

      std::vector<int> select;
      for (std::vector<StrawHitPtr>::const_iterator iTCHit=tclust._selectedTrackerHits.begin(); iTCHit!=tclust._selectedTrackerHits.end(); ++iTCHit){
	size_t index = iTCHit->key();
	select.push_back(index);
      }

      if (_doMinPrints) { std::cout<<"seed state "<<gRandom->GetSeed()<<std::endl; }
      std::sort(select.begin(),select.end());
      TObjArray *array=mu2eHits2ilc(shits,shitmcs,&select,-10,_useZCoordinate);
      //TObjArray *array=mu2eHits2ilc(shits,shitmcs,&select,tclust._meanTime-7.32505e+01-4.35,_useZCoordinate);
      itracker->LoadClusters(array);
      itracker->Clusters2Tracks();

      //fill seed information
      TrackSeed museed;
      museed._relatedTimeCluster=TrackerHitTimeClusterPtr(tclustHandle,ipeak);
      museed._t0    = tclust._minHitTime+6.59195;
      museed._errt0 = 5.;//tpeakerr;


      std::vector<int> usedhits;
      TObjArray *seeds=itracker->GetSeeds();

      std::vector<std::pair<double,int> > quals;
      for(int i=0;i<seeds->GetEntriesFast();i++){
	IlcDCHseed *pt=(IlcDCHseed*)seeds->At(i);
	if(!pt) continue;
	double qlt=pt->GetNumberOfClusters()-pt->GetChi2()/(pt->GetNumberOfClusters()-5);
	quals.push_back(std::pair<double,int>(qlt,i));
      }
      std::sort(quals.begin(),quals.end());
      bool seedfirst=true;
      for(int itr=quals.size()-1;itr>=0;itr--){
	IlcDCHseed *pt=(IlcDCHseed*)seeds->At(quals[itr].second);
	if(!pt) continue;
        if(pt->GetP()<0.05||pt->GetP()>0.15) continue;
        if(pt->Get1Pt()>0||pt->GetTgl()<0) continue;
        if(pt->GetNumberOfClusters()<15) continue;
        if(pt->GetChi2()/(pt->GetNumberOfClusters()-5)>10) continue;
	HelixVal segment;
	for(int ilay=pt->GetFirstPoint();ilay<pt->GetLastPoint();ilay++){
	  for(int j=0;j<kMaxInLayer;j++){
	    IlcDCHcluster*cl=pt->GetClusterPointer(ilay,j);
	    if(!cl)break;
	    usedhits.push_back(cl->GetIndex());
	    segment._selectedTrackerHitsIdx.push_back(cl->GetIndex());
	  }
	}
	double xy[2]={0,0};
	pt->PropagateToDCA(xy,IlcTracker::GetBz());
	pt->Print("params");
	segment._d0=pt->GetY()*10*pt->GetDir();
        segment._phi0=pt->GetAlpha()+(pt->GetDir()<0?TMath::Pi():0);
        segment._omega=pt->GetC()/10.;
        segment._z0=pt->GetZ()*10;
        segment._tanDip=pt->GetTgl();
	int pt2pt[5]={0,2,4,1,3};
	double coeff[5]={-10.*pt->GetDir(),1,0.1*IlcTracker::GetBz()*kB2C,10,1};
	for(int i=0;i<5;i++)
	  for(int j=0;j<5;j++)
	    segment._covMtrx[i][j]=pt->Cov(pt2pt[i],pt2pt[j])*coeff[i]*coeff[j];

	if(seedfirst){
	  museed._fullTrkSeed=segment;
          seedfirst=false;
	}
	museed._loopSeeds.push_back(segment);
      }
      std::sort(usedhits.begin(),usedhits.end());
      museed._fullTrkSeed._selectedTrackerHitsIdx.clear();

      int iprev=-1;
      int nunic=0;
      for(size_t i=0;i<usedhits.size();i++){
	if(usedhits[i]!=iprev) {
	  nunic++;
	  museed._fullTrkSeed._selectedTrackerHitsIdx.push_back( HitIndex(usedhits[i]) );
	  museed._selectedTrackerHits.push_back( StrawHitPtr (cellhits,usedhits[i]) );

	}
	iprev=usedhits[i];
      }
      if (_doMinPrints) { std::cout<<"used unic hits "<<nunic<<" sum "<<usedhits.size()<<" from "<<array->GetEntriesFast()<<std::endl; }

      if(nunic!=0)
	outseeds->push_back(museed);
      
      itracker->UnloadClusters();

      array->Delete();
      delete array;

      //      std::cout<<"------- start Tracks search for event "<<event.id().event()<<" ------"<<std::endl;


    }


    event.put(outseeds);

    clock_t stopClock = clock();
    if (_doMinPrints) { std::cout<<"-------- N clok to analyze 1 ev by ITTrackMiniIReco "<<stopClock-startClock<<" @ "<<CLOCKS_PER_SEC<<std::endl; }
    
  } // end produce


  void ITTrackMiniIReco::endJob(){
  }

}  // end namespace mu2e

using mu2e::ITTrackMiniIReco;
DEFINE_ART_MODULE(ITTrackMiniIReco);

