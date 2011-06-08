//
// Module to perform BaBar Kalman fit
//
// $Id: KalFitTest_module.cc,v 1.2 2011/06/08 03:53:11 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/08 03:53:11 $
//

// framework
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "BaBar/PdtPid.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/DetStrawHitElem.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "BField/BFieldFixed.hh"
#include "KalmanTrack/KalHit.hh"
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
// simple structs
  struct threevec {
    Float_t _x,_y,_z;
    threevec(): _x(0.0),_y(0.0),_z(0.0) {}
    threevec(const CLHEP::Hep3Vector& vec) : _x(vec.x()),_y(vec.y()),_z(vec.z()) {}
    threevec& operator = (const CLHEP::Hep3Vector& vec){ _x =vec.x(); _y =vec.y(); _z= vec.z(); return *this; }
  };
  
  struct helixpar {
    Float_t _d0, _p0, _om, _z0, _td;
    helixpar() : _d0(0.0),_p0(0.0),_om(0.0),_z0(0.0),_td(0.0) {}
    helixpar(const HepVector& pvec) : _d0(pvec[0]),_p0(pvec[1]),_om(pvec[2]),_z0(pvec[3]),_td(pvec[4]) {}
    helixpar(const HepSymMatrix& pcov) : _d0(sqrt(pcov.fast(1,1))),_p0(sqrt(pcov.fast(2,2))),_om(sqrt(pcov.fast(3,3))),
      _z0(sqrt(pcov.fast(4,4))),_td(sqrt(pcov.fast(5,5))) {}
  };
    
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
    // bfield object
    BField* _bfield;
    // configuration parameters
    int _diag;
    int _printfreq;
    std::string _strawhitslabel;
    // cache of event objects
    const StrawHitCollection* _strawhits;
    const StrawHitMCTruthCollection* _mcstrawhits;
    const PtrStepPointMCVectorCollection* _mchitptr;
    const StepPointMCCollection* _mcsteps;
    // fitter
    KalFit _kfit;
    // helper functions
    bool findData(art::Event& e);
    bool findMC(art::Event& e);
    bool trkFromMC(TrkDef& trkdef);
    void hitDiag(const TrkStrawHit* strawhit);
    void trkDiag(TrkDef const& mytrk, TrkKalFit const& myfit);
// trk tuple variables
    TTree *_trkdiag;
    UInt_t _eventid;
    UInt_t _trkid;
    Int_t _fitstatus;
    Float_t _t00;
    Float_t _t00err;
    Float_t _t0;
    Float_t _t0err;
    Float_t _mct0;
    Int_t _nhits;
    Int_t _ndof;
    UInt_t _niter;
    UInt_t _nt0iter;
    UInt_t _nweediter;
    Int_t _nactive;
    Float_t _chisq;
    Float_t _fitmom;
    Float_t _fitmomerr;
    Float_t _mcmom;
    Float_t _seedmom;
    
    helixpar _fitpar;
    helixpar _fiterr;
    helixpar _mcpar;    
    
// stl not supported by ancient versions of root: FIXME!!!!
//    std::vector<threevec> _shposs;

// hit tuple variables
    TTree *_hitdiag;
    threevec _shpos;
    Float_t _dmid;
    Float_t _dmiderr;
    Float_t _hitt0;
    Float_t _hitt0err;
    Float_t _rdrift;
    Float_t _rdrifterr;
    Float_t _resid;
    Float_t _residerr;
    Float_t _edep;
    Int_t _amb;
    Float_t _hflt;
    Float_t _trkflt;
    Bool_t _active;
    Int_t _use;
    UInt_t _nmcsteps;
    threevec _mcpos;
    Float_t _mcdmid;
    Float_t _mchitt0;
    Float_t _mcrdrift;
// general
    static const double _vlight;
    static const double _vdrift; 
// 
  };
// convert speed of light from m/sec to mm/nsec
// actually, contrary to the documentation, I find this is already 
  const double KalFitTest::_vlight = CLHEP::c_light;
// drift velocity should come from a service FIXME!!!  
  const double KalFitTest::_vdrift = 0.05; // 50 um/nsec
  
  
  KalFitTest::KalFitTest(fhicl::ParameterSet const& pset) : 
    _diag(pset.get<int>("diagLevel",0)),
    _printfreq(pset.get<int>("printFrequency",10)),
    _strawhitslabel(pset.get<std::string>("strawHitsLabel")),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit"))
  {
  }

  KalFitTest::~KalFitTest(){}

  void KalFitTest::beginJob( )
  {
    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      
      _trkdiag=tfs->make<TTree>("trkdiag","trk diagnostics");
      _trkdiag->Branch("eventid",&_eventid,"eventid/i");
      _trkdiag->Branch("trkid",&_trkid,"trkid/i");
      _trkdiag->Branch("fitstatus",&_fitstatus,"fitstatus/I");
      _trkdiag->Branch("t00",&_t00,"t00/F");
      _trkdiag->Branch("t00err",&_t0err,"t00err/F");
      _trkdiag->Branch("t0",&_t0,"t0/F");
      _trkdiag->Branch("t0err",&_t0err,"t0err/F");
      _trkdiag->Branch("mct0",&_mct0,"mct0/F");
      _trkdiag->Branch("nhits",&_nhits,"nhits/I");
      _trkdiag->Branch("ndof",&_ndof,"ndof/I");
      _trkdiag->Branch("niter",&_niter,"niter/i");
      _trkdiag->Branch("nt0iter",&_nt0iter,"nt0iter/i");
      _trkdiag->Branch("nweediter",&_nweediter,"nweediter/i");
      _trkdiag->Branch("nactive",&_nactive,"nactive/I");
      _trkdiag->Branch("chisq",&_chisq,"chisq/F");
      _trkdiag->Branch("fitmom",&_fitmom,"fitmom/F");
      _trkdiag->Branch("fitmomerr",&_fitmomerr,"fitmomerr/F");
      _trkdiag->Branch("seedmom",&_seedmom,"seedmom/F");
      _trkdiag->Branch("mcmom",&_mcmom,"mcmom/F");
      _trkdiag->Branch("fitpar",&_fitpar,"d0/F:p0/F:om/F:z0/F:td/F");
      _trkdiag->Branch("fiterr",&_fiterr,"d0err/F:p0err/F:omerr/F:z0err/F:tderr/F");
      _trkdiag->Branch("mcpar",&_mcpar,"mcd0/F:mcp0/F:mcom/F:mcz0/F:mctd/F");
      
//      _trkdiag->Branch("shposs",&_shposs);
      
      _hitdiag=tfs->make<TTree>("hitdiag","hit diagnostics");      
      _hitdiag->Branch("shpos",&_shpos,"x/F:y/F:z/F");
      _hitdiag->Branch("rdrift",&_rdrift,"rdrift/F");
      _hitdiag->Branch("rdrifterr",&_rdrifterr,"rdrifterr/F");
      _hitdiag->Branch("amb",&_amb,"amb/I");
      _hitdiag->Branch("hitt0",&_hitt0,"hitt0/F");
      _hitdiag->Branch("hitt0err",&_hitt0err,"hitt0err/F");
      _hitdiag->Branch("dmid",&_dmid,"dmid/F");
      _hitdiag->Branch("dmiderr",&_dmiderr,"dmiderr/F");
      _hitdiag->Branch("resid",&_resid,"resid/F");
      _hitdiag->Branch("residerr",&_residerr,"residerr/F");
      _hitdiag->Branch("hflt",&_hflt,"hflt/F");
      _hitdiag->Branch("trkflt",&_trkflt,"trkflt/F");
      _hitdiag->Branch("active",&_active,"active/B");
      _hitdiag->Branch("use",&_use,"use/I");
      _hitdiag->Branch("edep",&_edep,"edep/F");
      
      _hitdiag->Branch("nmcsteps",&_nmcsteps,"nmcsteps/i");
      _hitdiag->Branch("mcpos",&_mcpos,"x/F:y/F:z/F");
      _hitdiag->Branch("mcrdrift",&_mcrdrift,"mcrdrift/F");
      _hitdiag->Branch("mchitt0",&_mchitt0,"mchitt0/F");
      _hitdiag->Branch("mcdmid",&_mcdmid,"mcdmid/F");
      
//      _hitdiag->Branch("shdir","CLHEP::HepVector",&_shdir);
    }
  }

  void KalFitTest::beginRun(art::Run& ){

// get magnetic field map.  For now, fix this to a uniform field, should come from a
// field map service, FIXME!!!!
    _bfield=new BFieldFixed(0.,0.,1.,0.);
    _kfit.setField(_bfield);
  }

  void KalFitTest::produce(art::Event& event ) 
  {
// event printout
    int iev=event.id().event();
    if((iev%_printfreq)==0)cout<<"KalFitTest: event="<<iev<<endl;
    _eventid = iev;
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
// must initialize t0, momentum, initial trajectory.  These should come from patrec
// that doesn't yet exist. For now, take from the MC truth.  This also won't work when there's background, FIXME!!!!!
    TrkDef mytrk;
    if(trkFromMC(mytrk)){
// use this to create a track
      TrkKalFit myfit;
      _kfit.makeTrack(mytrk,myfit);
// test if fit succeeded
      if(myfit._fit.success()){
//  diagnostics
        if(_diag > 0){
          trkDiag(mytrk,myfit);
          if(_diag > 1){
            for(std::vector<TrkStrawHit*>::iterator ihit=myfit._hits.begin();ihit!=myfit._hits.end();ihit++){
              TrkStrawHit* trkhit = *ihit;
              hitDiag(trkhit);
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
    _mcstrawhits = 0; _mchitptr = 0; _mcsteps = 0;
    art::Handle<StrawHitMCTruthCollection> truthHandle;
    if(evt.getByLabel("makeSH",truthHandle))
      _mcstrawhits = truthHandle.product();
  // Get the persistent data about pointers to StepPointMCs
    art::Handle<PtrStepPointMCVectorCollection> _mchitptrHandle;
    if(evt.getByLabel("makeSH","StrawHitMCPtr",_mchitptrHandle))
      _mchitptr = _mchitptrHandle.product();
  // Get the persistent data about the StepPointMCs
    art::Handle<StepPointMCCollection> _mcstepsHandle;
    if(evt.getByLabel("g4run","tracker",_mcstepsHandle))
      _mcsteps = _mcstepsHandle.product();
    return _mcstrawhits != 0 && _mchitptr != 0 && _mcsteps != 0;
  }
  
// define seed helix and t0 using MC truth
  bool KalFitTest::trkFromMC(TrkDef& mytrk) {
// preset to failure
    bool retval(false);
// t0 is defined as the time when the particle passes through z=0.  Find the step points
// closest to this and interpolate.  This should be replaced by referencing the intersection
// with the null element at z=0, FIXME!!!!!
    double zmin(-1e6);
    double zmax(1e6);
    double tmin,tmax;
    unsigned nstraws = _strawhits->size();
    std::vector<size_t> indices;
    double mct0(0.0);
    CLHEP::Hep3Vector pos;
    CLHEP::Hep3Vector mom;
    for(unsigned istraw=0;istraw<nstraws;istraw++){
      PtrStepPointMCVector const& mcptr(_mchitptr->at(istraw));
      unsigned nprimary(0);
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = *mcptr[j];
  // make sure this is the primary particle: not sure if this works with bkgs merged FIXME!!!!
        if( mchit.trackId().asInt() != 1 ) continue;
        nprimary++;
        double z = mchit.position().z();
        if(z > 0 && z < zmax){
          zmax = z;
          tmax = mchit.time();
  // this t0 is a calorimeter artifact introduced in the tracking code, it must be subtracted.
  // eventually this part of the simulation should be moved to the calorimeter code, FIXME!!!
          mct0 = _mcstrawhits->at(istraw).t0();
        } else if (z < 0 && z > zmin){
          zmin = z;
          tmin = mchit.time();
          pos = mchit.position();
          mom = mchit.momentum();
        }
      }
      if(nprimary > 0)indices.push_back(istraw);
    }
    if(zmin > -1e5 && zmax < 1e5 ){
// interpolate linearly
// must correct for straw hit t0; not sure what that means physically????
      double t0 = tmin + fabs(zmin)*(tmax-tmin)/(zmax-zmin) + mct0;
// should set t0err according to average, but that's not working.  For now, set to fixed value.
// set to a negative value, to force t0 finding from data
      double t0err = -5.0;
// require a minimum momentum; should be a parameter, FIXME!!!
      static double minmom(60.0);
      if(mom.mag() > minmom){
        double charge(-1.); // not sure how to get this from MC: FIXME!!!!!
        HepVector parvec(5,0);
// Babar interface still uses HepPoint: FIXME!!!!
        double hflt(0.0);
        TrkHelixUtils::helixFromMom( parvec, hflt, 
          HepPoint(pos.x(),pos.y(),pos.z()),
          mom,charge,*_bfield);
  // dummy covariance matrix; this should be set according to measured values, FIXME!!!!!
        HepSymMatrix dummy(5,1); 
        dummy(1,1)=1.; dummy(2,2)=0.1*0.1;dummy(3,3)=1e-2*1e-2;
        dummy(4,4)=1.; dummy(5,5)=0.1*0.1;
        mytrk = TrkDef(_strawhits,indices,parvec,dummy,t0,t0err);
        retval = true;
      }
    }
    return retval;
  }
 
  void KalFitTest::hitDiag(const TrkStrawHit* strawhit) {
// straw info
    _edep = strawhit->strawHit().energyDep();
//    CLHEP::Hep3Vector dmid = strawhit->wirePosition()-strawhit->straw().getMidPoint(); 
// trkhit  info
    _rdrift = strawhit->driftRadius();
    _rdrifterr = strawhit->driftRadiusErr();
    _amb = strawhit->ambig();
    CLHEP::Hep3Vector shpos;
    strawhit->hitPosition(shpos);
    _shpos = shpos;
//    _shposs.push_back(shpos);
    _dmid = strawhit->timeDiffDist();
    _dmiderr = strawhit->timeDiffDistErr();
// residual info
    double resid,residerr;
    if(strawhit->resid(resid,residerr,true)){
      _resid = resid;
      _residerr = residerr;
    } else {
      _resid = _residerr = -100.0;
    }
    _hitt0 = strawhit->hitT0();
    _hitt0err = strawhit->hitT0Err();
    _hflt = strawhit->hitLen();
    _trkflt = strawhit->fltLen();
    _active = strawhit->isActive();
    _use = strawhit->usability();
// mc information
    unsigned istraw = strawhit->index();
    StrawHitMCTruth const& mcstrawhit = (_mcstrawhits->at(istraw));
    _mcrdrift = mcstrawhit.driftDistance();
    _mcdmid = mcstrawhit.distanceToMid();
// straw hit t0 has very strange properties: skip it!
//    _mchitt0 = mcstrawhit.t0();
    CLHEP::Hep3Vector mcpos;
    PtrStepPointMCVector const& mcptr(_mchitptr->at(istraw));
    _nmcsteps = mcptr.size();
    double esum(0.0);
    double mct0(0.0);
    for( size_t imc=0; imc< _nmcsteps; ++imc ) {
      StepPointMC const& mchit = *mcptr[imc];
      if( mchit.trackId().asInt() == 1 ){
        double edep = mchit.eDep();
        esum += edep;
        mcpos += mchit.position()*edep;
        mct0 += mchit.time()*edep;
      }
    }
    if(esum > 0.0){
      mcpos /= esum;
      mct0 /= esum;
    }
// must correct for straw hit t0; what does that represent physically????
    mct0 += mcstrawhit.t0();
    _mcpos = mcpos;
    _mchitt0 = mct0;
    _hitdiag->Fill();
  }
  
  void
  KalFitTest::trkDiag(TrkDef const& mytrk, TrkKalFit const& myfit) {
// mc t0 value
    _mct0 = mytrk.trkT0().t0();
// initial t0 value
    _t00 = myfit._t00.t0();
    _t00err = myfit._t00.t0Err();
// iterations
    _nt0iter = myfit._nt0iter;
    _nweediter = myfit._nweediter;
// final t0 value
    _t0 = myfit._t0.t0();
    _t0err = myfit._t0.t0Err();
    if(myfit._krep != 0 && myfit._krep->fitCurrent()){
      _trkid = myfit._krep->parentTrack()->id();
      _nhits = myfit._krep->hotList()->nHit();
      _fitstatus = myfit._krep->fitCurrent();
      _niter = myfit._krep->iterations();
      _ndof = myfit._krep->nDof();
      _nactive = myfit._krep->nActive();
      _chisq = myfit._krep->chisq();
      // get the fit at the first hit
      const TrkStrawHit* firsthit = dynamic_cast<const TrkStrawHit*>(myfit._krep->firstHit()->kalHit()->hitOnTrack());
      double fltlen=firsthit->fltLen();
      double loclen;
      const TrkSimpTraj* ltraj = myfit._krep->localTrajectory(fltlen,loclen);
      _fitpar = helixpar(ltraj->parameters()->parameter());
      _fiterr = helixpar(ltraj->parameters()->covariance());
      CLHEP::Hep3Vector fitmom = myfit._krep->momentum(fltlen);
      BbrVectorErr momerr = myfit._krep->momentumErr(fltlen);
      _fitmom = fitmom.mag();
      Hep3Vector momdir = fitmom.unit();
      HepVector momvec(3);
      for(int icor=0;icor<3;icor++)
        momvec[icor] = momdir[icor];
      _fitmomerr = sqrt(momerr.covMatrix().similarity(momvec));
      
      // get MC truth at this point
      unsigned istraw = firsthit->index();
      PtrStepPointMCVector const& mcptr(_mchitptr->at(istraw));
      _mcmom = -100;
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = *mcptr[j];
        if( mchit.trackId().asInt() == 1 ) {
          CLHEP::Hep3Vector mcmom = mchit.momentum();
          CLHEP::Hep3Vector mcpos = mchit.position();
          double mclen = loclen;
          HepVector mcpar(5);
          TrkHelixUtils::helixFromMom( mcpar, mclen, 
            HepPoint(mcpos.x(),mcpos.y(),mcpos.z()),
            mcmom,-1.,*_bfield);
          _mcmom = mchit.momentum().mag();
          _mcpar = helixpar(mcpar);
          break;
        }
      }
      CLHEP::Hep3Vector seedmom = TrkMomCalculator::vecMom(*(myfit._krep->seed()),*_bfield,0.0);
      _seedmom = seedmom.mag();      
    } else {
      _fitstatus = -1;
      _nhits = -1;
      _fitmom = -1.0;
    }
    _trkdiag->Fill(); 
  }
}

using mu2e::KalFitTest;
DEFINE_ART_MODULE(KalFitTest);
