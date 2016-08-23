//
// TTracker Kalman Fit launcher
//
// $Id: TrkRecFit_module.cc,v 1.4 2014/08/30 12:19:38 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/30 12:19:38 $
//
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"
// utiliites
#include "GeneralUtilities/inc/Angles.hh"
// data
#include "DataProducts/inc/Helicity.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
// Mu2e BaBar
#include "BTrkData/inc/TrkStrawHit.hh"
#include "TrkReco/inc/KalFit.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
// root 
#include "TH1F.h"
#include "TTree.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <float.h>
#include <vector>
using namespace std; 
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace mu2e 
{
  class TrkRecFit : public art::EDProducer
  {
    public:
      explicit TrkRecFit(fhicl::ParameterSet const&);
      virtual ~TrkRecFit();
      virtual void beginJob();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event ); 
      void endJob();
    private:
      unsigned _iev;
      // configuration parameters
      int _debug;
      int _printfreq;
      bool _addhits; 
      art::Handle<mu2e::StrawHitCollection> _strawhitsH;
      art::Handle<HelixSeedCollection> _trksSeedH;
      // event object labels
      string _shLabel;
      string _shfLabel;
      string _seedFinderLabel; // Label to the Pattern Recognition module.
      StrawHitFlag _addsel;
      StrawHitFlag _addbkg;
      double _maxdtmiss;
      // outlier cuts
      double _maxseeddoca,_maxhelixdoca,_maxadddoca, _maxaddchi;
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      Helicity _helicity; // cached value of helicity expected for this fit
      double _amsign; // cached value of angular momentum sign
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      StrawHitFlagCollection* _flags;
      const HelixSeedCollection * _hscol;
      // Kalman fitters.  Seed fit has a special configuration
      KalFit _seedfit, _kfit;
      PayloadSaver _payloadSaver;
      // helper functions
      bool findData(const art::Event& e);
      void filterOutliers(TrkDef& trkdef,double maxdoca);
      void findMissingHits(KalRep* krep, vector<hitIndex>& indices);
      bool Helix2Traj (const RobustHelix &helix, HelixTraj &traj);

      // flow diagnostic
      TH1F* _cutflow;
  };

  TrkRecFit::TrkRecFit(fhicl::ParameterSet const& pset) :
    _debug(pset.get<int>("debugLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _addhits(pset.get<bool>("addhits",true)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shfLabel(pset.get<string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _seedFinderLabel(pset.get<string>("SeedCollectionLabel","RobustHelixFinder")),
    _addsel(pset.get<vector<string> >("AddHitSelectionBits",vector<string>{} )),
    _addbkg(pset.get<vector<string> >("AddHitBackgroundBits",vector<string>{})),
    _maxdtmiss(pset.get<double>("DtMaxMiss",40.0)),
    _maxseeddoca(pset.get<double>("MaxSeedDoca",10.0)),
    _maxhelixdoca(pset.get<double>("MaxHelixDoca",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",4.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _seedfit(pset.get<fhicl::ParameterSet>("SeedFit",fhicl::ParameterSet()),_fdir),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet()),_fdir),
    _payloadSaver(pset)
  {
    produces<KalRepCollection>();
    produces<KalRepPtrCollection>();
    produces<KalRepPayloadCollection>();
    produces<StrawHitFlagCollection>();
  }

  TrkRecFit::~TrkRecFit(){}

  void TrkRecFit::beginJob(){
    // create a histogram of throughput: this is a basic diagnostic that should ALWAYS be on
    art::ServiceHandle<art::TFileService> tfs;
    _cutflow=tfs->make<TH1F>("cutflow","Cutflow",10,-0.5,9.5);
    _cutflow->GetXaxis()->SetBinLabel(1,"All Events");
    _cutflow->GetXaxis()->SetBinLabel(2,"Time Peak");
    _cutflow->GetXaxis()->SetBinLabel(3,"Helix Fit");
    _cutflow->GetXaxis()->SetBinLabel(4,"Seed Fit");
    _cutflow->GetXaxis()->SetBinLabel(5,"Kalman Fit");
  }

  void TrkRecFit::beginRun(art::Run& ){
 // calculate the helicity
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    // helicity is a purely geometric quantity, however it's easiest
    // to determine it from the kinematics (angular momentum and Z momentum)
    _amsign = copysign(1.0,-_tpart.charge()*field.z());
    _helicity = Helicity(static_cast<float>(_fdir.dzdt()*_amsign));
  }

  void TrkRecFit::produce(art::Event& event ) {
    _cutflow->Fill(0.0);
    // create output
    unique_ptr<KalRepCollection>    tracks(new KalRepCollection );
    unique_ptr<KalRepPtrCollection> trackPtrs(new KalRepPtrCollection );
    art::ProductID kalRepsID(getProductID<KalRepCollection>(event));
    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"TrkRecFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TrkRecFit: data missing or incomplete"<< endl;
    }
    // copy in the existing flags
    _flags = new StrawHitFlagCollection(*_shfcol);
    unique_ptr<StrawHitFlagCollection> flags(_flags );

    // loop over the Helices
    if(_hscol->size()>0)_cutflow->Fill(1.0);

    for (size_t iseed=0; iseed<_hscol->size(); ++iseed) {
      HelixSeed const& hseed(_hscol->at(iseed));
      // convert the helix to a fit trajectory.  This accounts for the physical particle direction
      HelixTraj htraj(TrkParams(HelixTraj::NHLXPRM));
      if(Helix2Traj(hseed._helix,htraj)){
	_cutflow->Fill(2.0);
  // create the track definition
	TrkDef seeddef(hseed._timeCluster,htraj,_tpart,_fdir);
  // filter outliers; this doesn't use drift information, just straw positions
	filterOutliers(seeddef,_maxhelixdoca);
      // now, fit the seed helix from the filtered hits
	KalRep *seedrep(0);
	_seedfit.makeTrack(_shcol,seeddef,seedrep);
	if(seedrep != 0 && seedrep->fitStatus().success()){
	  _cutflow->Fill(3.0);
        // extract the helix trajectory from the helix fit, and initialize the full Kalman fit with this
	  double locflt;
	  const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(seedrep->localTrajectory(seedrep->flt0(),locflt));
	  TrkDef kaldef(seeddef);
	  kaldef.setHelix(*htraj);
        // filter the outliers
	  filterOutliers(kaldef,_maxseeddoca);
	  KalRep *krep(0);
	  _kfit.makeTrack(_shcol,kaldef,krep);
        // if successfull, try to add missing hits
	  if(krep != 0 && krep->fitStatus().success()){
	    _cutflow->Fill(4.0);
	    if(_addhits){
            // first, add back the hits on this track
	      _kfit.unweedHits(krep,_maxaddchi);
	      vector<hitIndex> misshits;
	      findMissingHits(krep,misshits);
	      if(misshits.size() > 0){
		_kfit.addHits(krep,_shcol,misshits,_maxaddchi);
	      }
	    }
	    if(krep != 0 && krep->fitStatus().success()){
	      if(!krep->fitCurrent()){
		cout << "Fit not current! " << endl;
	      }
	  // flag the hits used in this track. 
	      if(iseed<StrawHitFlag::_maxTrkId){
		for(auto ihit=krep->hitVector().begin();ihit != krep->hitVector().end();++ihit){
		  if((*ihit)->isActive())_flags->at(static_cast<TrkStrawHit*>(*ihit)->index()).merge(StrawHitFlag::trackBit(iseed));
		}
	      }
	// save successful kalman fits in the event
	      tracks->push_back(krep);
	      int index = tracks->size()-1;
	      trackPtrs->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
	    } else
	      delete krep;
	  }
	}
      // cleanup the seed fit
	delete seedrep;
      }
    }
    // put the tracks into the event
    art::ProductID tracksID(getProductID<KalRepPayloadCollection>(event));
    _payloadSaver.put(*tracks, tracksID, event);
    event.put(move(tracks));
    event.put(move(trackPtrs));
    event.put(move(flags));
  }

  void TrkRecFit::endJob(){
    // does this cause the file to close?
    art::ServiceHandle<art::TFileService> tfs;
  }

  // find the input data objects 
  bool TrkRecFit::findData(const art::Event& evt){
    _shcol = 0;
    _shfcol = 0;
    _hscol = 0;

    if(evt.getByLabel(_shLabel,_strawhitsH))
      _shcol = _strawhitsH.product();
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(evt.getByLabel(_shfLabel,shflagH))
      _shfcol = shflagH.product();
    if (evt.getByLabel(_seedFinderLabel,_trksSeedH))
      _hscol = _trksSeedH.product();
    return _shcol != 0 && _shfcol != 0 && _hscol != 0;
  }

  void TrkRecFit::filterOutliers(TrkDef& mydef,double maxdoca){
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    mydef.helix().getInfo(0.0,tpos,tdir);
    // tracker and conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    const vector<hitIndex>& indices = mydef.strawHitIndices();
    vector<hitIndex> goodhits;
    for(unsigned ihit=0;ihit<indices.size();++ihit){
      StrawHit const& sh = _shcol->at(indices[ihit]._index);
      Straw const& straw = tracker.getStraw(sh.strawIndex());
      CLHEP::Hep3Vector hpos = straw.getMidPoint();
      CLHEP::Hep3Vector hdir = straw.getDirection();
      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
      TrkLineTraj htraj(spt,hdir,-20,20);
      // estimate flightlength along track.  This assumes a constant BField!!!
      double fltlen = (hpos.z()-tpos.z())/tdir.z();
      TrkPoca hitpoca(mydef.helix(),fltlen,htraj,0.0);
      // keep hits with small residuals
      if(fabs(hitpoca.doca()) < maxdoca){
        goodhits.push_back(indices[ihit]);
      }
    }
    // update track
    mydef.setIndices(goodhits);
  }

  void TrkRecFit::findMissingHits(KalRep* krep,vector<hitIndex>& misshits) {
    const Tracker& tracker = getTrackerOrThrow();
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tpos;
    krep->pieceTraj().getInfo(0.0,tpos,tdir);
    unsigned nstrs = _shcol->size();
    TrkStrawHitVector tshv;
    convert(krep->hitVector(),tshv);
    for(unsigned istr=0; istr<nstrs;++istr){
      if(_flags->at(istr).hasAllProperties(_addsel)&& !_flags->at(istr).hasAnyProperty(_addbkg)){
        StrawHit const& sh = _shcol->at(istr);
        if(fabs(_shcol->at(istr).time()-krep->t0()._t0) < _maxdtmiss) {
          // make sure we haven't already used this hit

          vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
          if(ifnd == tshv.end()){
            // good in-time hit.  Compute DOCA of the wire to the trajectory
            Straw const& straw = tracker.getStraw(sh.strawIndex());
            CLHEP::Hep3Vector hpos = straw.getMidPoint();
            CLHEP::Hep3Vector hdir = straw.getDirection();
            // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
            HepPoint spt(hpos.x(),hpos.y(),hpos.z());
            TrkLineTraj htraj(spt,hdir,-20,20);
            // estimate flightlength along track.  This assumes a constant BField!!!
            double fltlen = (hpos.z()-tpos.z())/tdir.z();
            TrkPoca hitpoca(krep->pieceTraj(),fltlen,htraj,0.0);
            // flag hits with small residuals
            if(fabs(hitpoca.doca()) < _maxadddoca){
              misshits.push_back(istr);
            }
          }
        }
      }
    }
  }

  bool TrkRecFit::Helix2Traj (RobustHelix const& helix, HelixTraj &traj) {
    // compare the input with this configuration's helicity: these must be the same
    bool retval = _helicity == helix.helicity();
    // radius = 0 or lambda=0 are degenerate cases that this representation can't handle
    retval &= helix.radius() > 0.0 && helix.lambda() != 0.0;
    if(retval){
      HepVector helParams(5);
      // radius and omega have inverse magnitude, omega is signed by the angular momentum 
      helParams[HelixTraj::omegaIndex] = _amsign/helix.radius();
      // phi0 is the azimuthal angle of the particle velocity vector at the point
      // of closest approach to the origin.  It's sign also depends on the angular
      // momentum.  To translate from the center, we need to reverse coordinates
      helParams[HelixTraj::phi0Index] = atan2(-_amsign*helix.center().x(),_amsign*helix.center().y());
      // d0 describes the distance to the origin at closest approach.
      // It is signed by the particle angular momentum WRT the origin.
      // The Helix fit radial bias is anti-correlated with d0; correct for it here.
      helParams[HelixTraj::d0Index] = _amsign*(helix.center().perp() - helix.radius());
      // the dip angle is measured WRT the perpendicular, signed by the z component of linear momentum
      helParams[HelixTraj::tanDipIndex] = _amsign*helix.lambda()/helix.radius();
      // must change conventions here: fz0 is the phi at z=0, z0 is defined at the point of closest approach
      // resolve the loop ambiguity such that the POCA is closest to z=0.
      double phi = helix.fz0()+_amsign*M_PI_2;
      double dphi = Angles::deltaPhi(phi,helParams[HelixTraj::phi0Index]);
      // choose z0 (which loop) so that f=0 is as close to z=0 as possible
      helParams[HelixTraj::z0Index] = dphi*helParams[HelixTraj::tanDipIndex]/helParams[HelixTraj::omegaIndex]; 
      // setup a dummy error matrix.
      HepVector perr(5,0);
      // estimated parameter errors based on average performance.  These should be parameters, FIXME!!!
      // These get scaled up by the fit in early iterations
      perr[HelixTraj::d0Index]     = 34.0;
      perr[HelixTraj::phi0Index]   = 0.02;
      perr[HelixTraj::omegaIndex]  = 0.0002;
      perr[HelixTraj::tanDipIndex] = 0.05;
      perr[HelixTraj::z0Index]     = 15.0;

      CLHEP::HepSymMatrix covar = vT_times_v(perr);
      traj = HelixTraj(helParams,covar);
    }
    return retval;
  }

}
using mu2e::TrkRecFit;
DEFINE_ART_MODULE(TrkRecFit);
