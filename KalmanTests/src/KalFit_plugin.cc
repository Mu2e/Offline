//
// Module to perform BaBar Kalman fit
//
// $Id: KalFit_plugin.cc,v 1.1 2011/06/02 00:00:05 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/06/02 00:00:05 $
//

// framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
// data
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHit.hh"
#include "ToyDP/inc/StrawHitMCTruth.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "BaBar/PdtPid.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/DetStrawHitElem.hh"

#include "KalmanTrack/KalContext.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalHit.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkPoca.hh"
#include "BaBar/ErrLog.hh"
#include "BField/BField.hh"
#include "BField/BFieldFixed.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetMaterial.hh"
#include "MatEnv/MatDBInfo.hh"
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
#include "FWCore/Framework/interface/EDProducer.h"


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
  
  struct fltlencomp : public binary_function<TrkStrawHit*, TrkStrawHit*, bool> {
	  bool operator()(TrkStrawHit* x, TrkStrawHit* y) { return x->fltLen() < y->fltLen(); }
  };
  
  class KalFit : public edm::EDProducer
  {
  public:
    explicit KalFit(edm::ParameterSet const&);
    virtual ~KalFit();
    virtual void beginJob(edm::EventSetup const&);
    virtual void beginRun(edm::Run&,edm::EventSetup const&);
    void endJob();
    void produce(edm::Event& e, edm::EventSetup const&);

  private:

    KalContext* _kalcon; // BaBar configuration object
    BField* _bfield; // magentic field description, BaBar wrapper class
    DetStrawHitElem _wallelem; // fake element to represent straw hit material
    DetStrawHitElem _gaselem; // fake element to represent straw hit material
    // configuration parameters
    int _debug;
    int _diag;
    bool _fieldcorr;
    bool _material;
    bool _ambigflip;
    bool _weedhits;
    bool _updatet0;
    double _t0tol;
    double _maxhitchi;
    unsigned _maxiter;
    unsigned _printfreq;
    unsigned _minnstraws;
    unsigned _maxweed;
    std::string _strawhitslabel;
    // cache of event objects
    const StrawHitCollection* _strawhits;
    const StrawHitMCTruthCollection* _mcstrawhits;
    const DPIndexVectorCollection* _mchitptr;
    const StepPointMCCollection* _mcsteps;
    // helper functions
    bool findData(edm::Event& e);
    bool findMC(edm::Event& e);
    bool fitable(const TrkDef& mytrk);
    bool updateT0(KalRep* kalrep);
    TrkErrCode weedHits(KalRep* kalrep, unsigned& niter);
    bool trkFromMC(TrkDef& trkdef);
    void makeHits(TrkDef& mytrk, std::vector<TrkStrawHit*>& hits);
    void hitDiag(const TrkStrawHit* strawhit);
    void trkDiag(const KalRep* kalrep);
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
  const double KalFit::_vlight = CLHEP::c_light;
// drift velocity should come from a service FIXME!!!  
  const double KalFit::_vdrift = 0.05; // 50 um/nsec
  
// ugly place to put this
  MatDBInfo* DetStrawHitType::_matdbinfo(new MatDBInfo);  
  
  KalFit::KalFit(edm::ParameterSet const& pset) : _wallelem(wall),_gaselem(gas),
    _debug(pset.getUntrackedParameter<int>("debugLevel",0)),
    _diag(pset.getUntrackedParameter<int>("diagLevel",0)),
    _fieldcorr(pset.getUntrackedParameter<bool>("fieldCorrections",false)),
    _material(pset.getUntrackedParameter<bool>("material",false)),
    _ambigflip(pset.getUntrackedParameter<bool>("ambigflip",false)),
    _weedhits(pset.getUntrackedParameter<bool>("weedhits",true)),
    _updatet0(pset.getUntrackedParameter<bool>("updateT0",true)),
    _t0tol(pset.getUntrackedParameter<double>("t0Tolerance",1.0)),
    _maxhitchi(pset.getUntrackedParameter<double>("maxhitchi",5.0)),
    _maxiter(pset.getUntrackedParameter<unsigned>("maxiter",3)),
    _printfreq(pset.getUntrackedParameter<unsigned>("printfrequency",1)),
    _minnstraws(pset.getUntrackedParameter<unsigned>("minnstraws",20)),
    _maxweed(pset.getUntrackedParameter<unsigned>("maxweed",10)),
    _strawhitslabel(pset.getParameter<std::string>("strawHitsLabel"))
  {
  }

  KalFit::~KalFit(){}

  void KalFit::beginJob(edm::EventSetup const& )
  {
    if(_diag > 0){
      edm::Service<edm::TFileService> tfs;
      
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

  void KalFit::beginRun(edm::Run&,edm::EventSetup const& ){    
// get magnetic field map.  For now, fix this to a uniform field, should come from a
// field map service, FIXME!!!!
    _bfield=new BFieldFixed(0.,0.,1.,0.);
    double Bval = _bfield->bFieldNominal();
    if(_debug > 0)cout<<"BVal "<<Bval<<" GeV/c/Tesla/cm "<<BField::cmTeslaToGeVc<<endl;
    _kalcon = new KalContext;
    _kalcon->setBendSites(_fieldcorr);
    _kalcon->setMaterialSites(_material);
    _kalcon->setForbidAmbigFlips(_ambigflip); //false: free left-rigth ambiguity, true: will be fixed from sim
    _kalcon->setMaxIterations(_maxiter);
    // these are currently fixed, they should be set as parameters and re-optimized FIXME!!!!
    _kalcon->setMaxIntersections(0);
    _kalcon->setMaxDMom(10);
    _kalcon->setSmearFactor(1e6);
    _kalcon->setMinDOF(15,TrkEnums::bothView);
    _kalcon->setMinDOF(20,TrkEnums::xyView);
    _kalcon->setMinDOF(0,TrkEnums::zView);
    _kalcon->setIntersectionTolerance(100);
    _kalcon->setMaxMomDiff(1.0); // 1 MeV
    _kalcon->setTrajBuffer(0.01); // 10um
    _kalcon->setMinGap(0.0); // no minimum separation between sites
    _kalcon->setDefaultType(PdtPid::electron); // by default, fit electrons
    
  }

  void KalFit::produce(edm::Event& event, edm::EventSetup const&) 
  {
// event printout
    int iev=event.id().event();
    if((iev%_printfreq)==0)cout<<"KalFit: event="<<iev<<endl;
    _eventid = iev;
// find the data
    if(!findData(event)){
      cout << "No straw hits found " << endl;
      return;
    }
// find mc truth
    if(_diag > 0){
      if(!findMC(event)){
        cout<<"MC information missing "<< endl;
// for now, this is a fatal error
        return;
      }
    }
// must initialize t0, momentum, initial trajectory.  These should come from patrec
// that doesn't yet exist. For now, take from the MC truth.  This also won't work when there's background, FIXME!!!!!
    TrkDef mytrk;
    if(trkFromMC(mytrk)){
// test if fitable
      if(fitable(mytrk)){
// remember mc t0 value
        _mct0 = mytrk.t0();
// create the track hits from the straw hits
        std::vector<TrkStrawHit*> hits;
        makeHits(mytrk,hits);
// remember initial t0 value
        _t00 = mytrk.t0();
        _t00err = mytrk.t0Err();
// Create the BaBar hit list.  This takes ownership
// Also create a straw hit intersection for each straw hit (active or not)
        std::vector<DetIntersection> detinter;
        TrkHotListFull* hotlist = new TrkHotListFull();
        for(std::vector<TrkStrawHit*>::iterator ihit=hits.begin();ihit!=hits.end();ihit++){
          TrkStrawHit* trkhit = *ihit;
          hotlist->append(trkhit);
          double fltlen = trkhit->fltLen();
// note this accounts for the material of both walls
          double wallpath = trkhit->wallPath();
          double gaspath = trkhit->gasPath();
  // offset the paths to avoid stacking elements
          double wlen = fltlen - 0.5*trkhit->straw().getRadius();
          double glen = fltlen + 0.5*trkhit->straw().getRadius();
          detinter.push_back(DetIntersection(&_wallelem,&mytrk.helix(),
            wlen,wlen-wallpath,wlen+wallpath));
          detinter.push_back(DetIntersection(&_gaselem,&mytrk.helix(),
            glen,glen-gaspath,glen+gaspath));
        }
// Create BaBar track and Kalman fit
        TrkRecoTrk* trk = new TrkRecoTrk(_kalcon->defaultType(), 0, 0);
        assert(trk != 0);
        trk->setBField(_bfield);
// create Kalman rep
        KalRep *kalrep = new KalRep(mytrk.helix(), hotlist, detinter, trk, *_kalcon, PdtPid::electron);
        assert(kalrep != 0);
        trk->setRep(kalrep);
        trk->resetT0(mytrk.t0());
// fit the track
        TrkErrCode err = kalrep->fit();
// update t0, and propagate it to the hits
        double oldt0 = 0;
        _nt0iter = 0;
        while(err.success() && _updatet0 && fabs(kalrep->trackT0()-oldt0) > _t0tol && 
          _nt0iter < _kalcon->maxIterations()){
          oldt0 = kalrep->trackT0();
          if(updateT0(kalrep)){
            kalrep->resetFit();
            err = kalrep->fit();
            _nt0iter++;
          } else
            break;
// drop outlyers
          if(_weedhits){
            _nweediter = 0;
            err = weedHits(kalrep,_nweediter);
          }
        }        
//  diagnostics
        if(_diag > 0){
          trkDiag(kalrep);
          for(std::vector<TrkStrawHit*>::iterator ihit=hits.begin();ihit!=hits.end();ihit++){
            TrkStrawHit* trkhit = *ihit;
            hitDiag(trkhit);
          }
        }
// cleanup
        delete trk;
      }
    }
  }
  
  void KalFit::endJob()
  {
    edm::Service<edm::TFileService> tfs;
  }
  
// find the input data objects needed to make tracks
  bool KalFit::findData(edm::Event& evt){
    _strawhits = 0;
    edm::Handle<mu2e::StrawHitCollection> strawhitsH;
    if(evt.getByLabel(_strawhitslabel,strawhitsH))
      _strawhits = strawhitsH.product();
    return _strawhits != 0;
  }
  
// find the MC truth objects in the event and set the local cache
  bool KalFit::findMC(edm::Event& evt) {
    _mcstrawhits = 0; _mchitptr = 0; _mcsteps = 0;
    edm::Handle<StrawHitMCTruthCollection> truthHandle;
    if(evt.getByLabel("makeSH",truthHandle))
      _mcstrawhits = truthHandle.product();
  // Get the persistent data about pointers to StepPointMCs
    edm::Handle<DPIndexVectorCollection> _mchitptrHandle;
    if(evt.getByLabel("makeSH","StrawHitMCPtr",_mchitptrHandle))
      _mchitptr = _mchitptrHandle.product();
  // Get the persistent data about the StepPointMCs
    edm::Handle<StepPointMCCollection> _mcstepsHandle;
    if(evt.getByLabel("g4run","tracker",_mcstepsHandle))
      _mcsteps = _mcstepsHandle.product();
    return _mcstrawhits != 0 && _mchitptr != 0 && _mcsteps != 0;
  }
  
  bool
  KalFit::fitable(const TrkDef& mytrk){
    return mytrk.strawHitIndices().size() >= _minnstraws;
  }
  
// define seed helix and t0 using MC truth
  bool KalFit::trkFromMC(TrkDef& mytrk) {
// preset to failure
    bool retval(false);
// t0 is defined as the time when the particle passes through z=0.  Find the step points
// closest to this and interpolate
    int imin(-1), imax(-1);
    double zmin(-1e6);
    double zmax(1e6);
    unsigned nstraws = _strawhits->size();
    std::vector<size_t> indices;
    double mct0(0.0);
    for(unsigned istraw=0;istraw<nstraws;istraw++){
      DPIndexVector const& mcptr(_mchitptr->at(istraw));
      unsigned nprimary(0);
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = (*_mcsteps)[mcptr[j].index];
  // make sure this is the primary particle: not sure if this works with bkgs merged FIXME!!!!
        if( mchit.trackId().asInt() != 1 ) continue;
        nprimary++;
        double z = mchit.position().z();
        if(z > 0 && z < zmax){
          zmax = z;
          imax = mcptr[j].index;
          mct0 = _mcstrawhits->at(istraw).t0();
        } else if (z < 0 && z > zmin){
          zmin = z;
          imin = mcptr[j].index;
        }
      }
      if(nprimary > 0)indices.push_back(istraw);
    }
    if(imin >= 0 && imax >= 0){
// interpolate linearly
      double tmin = (*_mcsteps)[imin].time();
      double tmax = (*_mcsteps)[imax].time();
// must correct for straw hit t0; not sure what that means physically????
      double t0 = tmin + fabs(zmin)*(tmax-tmin)/(zmax-zmin) + mct0;
// should set t0err according to average, but that's not working.  For now, set to fixed value.
// set to a negative value, to force t0 finding from data
      double t0err = -10.0;
// get the position and momentum at the point nearest the center to define the helix
      unsigned ibest = fabs(zmax) < fabs(zmin) ? imax : imin;
      CLHEP::Hep3Vector pos=(*_mcsteps)[ibest].position();
      CLHEP::Hep3Vector mom=(*_mcsteps)[ibest].momentum();
// require a minimum momentum; should be a parameter, FIXME!!!
      static double minmom(60.0);
      if(mom.mag() > minmom){
        double charge(-1.); // not sure how to get this from MC: FIXME!!!!!
        if(_debug > 1) {
          cout<<"pos= "<<pos<<", mom= "<<mom<<endl;
        }
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
  
  void
  KalFit::makeHits(TrkDef& mytrk, std::vector<TrkStrawHit*>& hits) {
    const Tracker& tracker = getTrackerOrThrow();
    // find flightlength at z=0
    double flt0 = mytrk.helix().zFlight(0.0);
    unsigned nind = mytrk.strawHitIndices().size();
    double tsum(0.0);
    for(unsigned iind=0;iind<nind;iind++){
      unsigned istraw = mytrk.strawHitIndices()[iind];
      const StrawHit& strawhit(mytrk.strawHitCollection()->at(istraw));
      const Straw& straw = tracker.getStraw(strawhit.strawIndex());
    // compute initial flightlength from helix and hit Z
      double hflt = mytrk.helix().zFlight(straw.getMidPoint().z());
    // estimate the time the track reaches this hit, assuming speed-of-light travel along the helix. Should use actual
    // speed based on momentum and assumed particle species FIXME!!!
      double tprop = (hflt -flt0)/_vlight;
      double hitt0 = mytrk.t0() + tprop;
    // subtract the propagation time and the average wire signal delay when computing hit time
    // using vlight = vwire, FIXME!!!!!
      tsum += strawhit.time() - tprop - straw.getHalfLength()/_vlight;
    // create the hit object
      TrkStrawHit* trkhit = new TrkStrawHit(strawhit,straw,istraw,hitt0,mytrk.t0Err());
      assert(trkhit != 0);
    // refine the flightlength, as otherwise hits in the same plane are at exactly the same flt, which can cause problems
      TrkPoca poca(mytrk.helix(),hflt,*trkhit->hitTraj(),trkhit->hitLen());
    // set the initial flightlength; this should be a constructor parameter, FIXME!!!
      if(poca.status().success()){
        trkhit->setFltLen(poca.flt1());
        trkhit->setHitLen(poca.flt2());
      }
      hits.push_back(trkhit);
    }
  // if the initial t0error was 0, compute t0 and override
    if(mytrk.t0Err() < 0.0 && nind > 0){
  // assuming a flat drift time means correcting by half the maximum drift time
      double tmax = hits[0]->straw().getRadius()/_vdrift;
      double t0 = tsum/nind - 0.5*tmax;
  // estimate the error using the same assumption
      double t0err = tmax/sqrt(12*nind);
  // now update the hits and trackdef object
      mytrk.setT0(t0,t0err);
      for(std::vector<TrkStrawHit*>::iterator ihit= hits.begin();ihit != hits.end(); ihit++){
        double hitt0 = t0 + ((*ihit)->fltLen() -flt0)/_vlight;
        (*ihit)->updateT0(hitt0,t0err);
      }
    }
  // sort the hits by flightlength
    std::sort(hits.begin(),hits.end(),fltlencomp());
  }
  
  void KalFit::hitDiag(const TrkStrawHit* strawhit) {
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
    DPIndexVector const& mcptr(_mchitptr->at(istraw));
    _nmcsteps = mcptr.size();
    double esum(0.0);
    double mct0(0.0);
    for( size_t imc=0; imc< _nmcsteps; ++imc ) {
      StepPointMC const& mchit = (*_mcsteps)[mcptr[imc].index];
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
  KalFit::trkDiag(const KalRep* kalrep) {
    _trkid = kalrep->parentTrack()->id();
    _nhits = kalrep->hotList()->nHit();
    _fitstatus = kalrep->fitCurrent();
    _niter = kalrep->iterations();
    _ndof = kalrep->nDof();
    _nactive = kalrep->nActive();
    _chisq = kalrep->chisq();
    _t0 = kalrep->trackT0();
    if(kalrep->fitCurrent()){
      // get the fit at the first hit
      const TrkStrawHit* firsthit = dynamic_cast<const TrkStrawHit*>(kalrep->firstHit()->kalHit()->hitOnTrack());
      double fltlen=firsthit->fltLen();
      double loclen;
      const TrkSimpTraj* ltraj = kalrep->localTrajectory(fltlen,loclen);
      _t0err = firsthit->hitT0Err();
      
      _fitpar = helixpar(ltraj->parameters()->parameter());
      _fiterr = helixpar(ltraj->parameters()->covariance());
            
      CLHEP::Hep3Vector fitmom = kalrep->momentum(fltlen);
      BbrVectorErr momerr = kalrep->momentumErr(fltlen);
      _fitmom	= fitmom.mag();
      Hep3Vector momdir = fitmom.unit();
      HepVector momvec(3);
      for(int icor=0;icor<3;icor++)
        momvec[icor] = momdir[icor];
      _fitmomerr = sqrt(momerr.covMatrix().similarity(momvec));
      
      // get MC truth at this point
      unsigned istraw = firsthit->index();
      DPIndexVector const& mcptr(_mchitptr->at(istraw));
      _mcmom = -100;
      for( size_t j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = (*_mcsteps)[mcptr[j].index];
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
      
    } else {
      _fitmom = -1.0;
    }
    CLHEP::Hep3Vector seedmom = TrkMomCalculator::vecMom(*(kalrep->seed()),*_bfield,0.0);
    _seedmom = seedmom.mag();
    _trkdiag->Fill();
    
  }
  
  bool
  KalFit::updateT0(KalRep* kalrep){
    bool retval(false);
// need to have a valid fit
    if(kalrep->fitValid()){
// update the fit so the hots are using the most recent trajectory
//      TrkErrCode update = kalrep->updateSites();
// find the global fltlen associated with z=0.  This should be a piectraj function, FIXME!!!
      double loclen;
      double flt0 = 0.0;
      double dz(10.0);
      unsigned niter(0);
      while(fabs(dz) > 1.0 && niter < _kalcon->maxIterations() ) {
        const HelixTraj* helix = dynamic_cast<const HelixTraj*>(kalrep->localTrajectory(flt0,loclen));
        flt0 += helix->zFlight(0.0)-loclen;
        dz = kalrep->traj().position(flt0).z();
        niter++;
      }
// find hits
      TrkHotList* hits = kalrep->hotList();
      std::vector<double> hitst0; // store t0, to allow outlyer removal
      for(TrkHotList::nc_hot_iterator ihit = hits->begin(); ihit != hits->end(); ihit++){
        TrkStrawHit* hit = dynamic_cast<TrkStrawHit*>(ihit.get());
        assert(hit != 0);
        if(hit->isActive() && hit->poca()!= 0 && hit->poca()->status().success()){
// copy the seed
          static TrkSimpTraj* straj = kalrep->seed()->clone();
// find the hit site in the rep
          const KalHit* hitsite = kalrep->findHotSite(hit);
// set helix to the local parameters EXCLUDING THIS HIT
          if(hitsite != 0 && kalrep->smoothedTraj(hitsite,straj)){
            TrkPoca poca(*straj,hit->fltLen(),*(hit->hitTraj()),hit->hitLen());
            if(poca.status().success()){
              double doca = fabs(poca.doca());
// require a minimum doca to avoid ambiguity bias.  mindoca shoudl be a parameter, FIXME!!!
              static double mindoca(0.4);
              if(doca > mindoca){
// propagation time to this hit from z=0.  This assumes beta=1, FIXME!!!
                double tflt = (hit->fltLen()-flt0)/_vlight;
// drift time of this hit (plus t0)
                double tdrift = hit->time() - tflt;
// t0 = Time difference between the drift time and the DOCA time.  sign of DOCA is irrelevant here.
                double hitt0 = tdrift - doca/_vdrift;
                hitst0.push_back(hitt0);
              }
            }
          }
        }
      }
      if(hitst0.size() >1){
// iterate over outlyer removal.  T0 window should be a parameter, FIXME!!!!
        double nsig(2.5);
        bool changed(true);
        double t0 = kalrep->trackT0();
        double t0err(-1.0);
        
        std::vector<bool> used(hitst0.size(),true);
        unsigned niter(0);
        while(changed && niter < 10){
          niter++;
          unsigned nactive(0);
          double t0sum(0.0);
          double t0sum2(0.0);
          for(unsigned ihit=0;ihit<hitst0.size();ihit++){
            if(used[ihit]){
              nactive++;
              t0sum += hitst0[ihit];
              t0sum2 += hitst0[ihit]*hitst0[ihit];
            }
          }
          t0 = t0sum/nactive;
          double t02 = t0sum2/nactive;
          double t0sig = sqrt(max(t02 - t0*t0,0.0));
          t0err = t0sig/sqrt(nactive);
          changed = false;
          for(unsigned ihit=0;ihit<hitst0.size();ihit++){
            bool useit = fabs(hitst0[ihit]-t0) < nsig*t0sig;
            changed |= useit != used[ihit];
            used[ihit] = useit;
          }
        }
          
// reset track t0; should also set t0err, FIXME!!!
        kalrep->parentTrack()->resetT0(t0);
// reset all the hit times
        for(TrkHotList::nc_hot_iterator ihit = hits->begin(); ihit != hits->end(); ihit++){
          TrkStrawHit* hit = dynamic_cast<TrkStrawHit*>(ihit.get());
          assert(hit != 0);
// correct for flightlength.  Again assumes beta=1, FIXME!!!
          double hitt0 = t0 + (hit->fltLen()-flt0)/_vlight;
          hit->updateT0(hitt0,t0err);
        }
        retval = true;
      }
    }
    return retval;
  }
  
  TrkErrCode
  KalFit::weedHits(KalRep* kalrep, unsigned& niter) {
    TrkErrCode retval;
    // Loop over HoTs and find HoT with largest contribution to chi2.  If this value
    // is greater than some cut value, deactivate that HoT and reFit

    double worst = -1.;
    TrkHitOnTrk* worstHot = 0;
    TrkHotList* hots = kalrep->hotList();
    for (TrkHotList::nc_hot_iterator iHot = hots->begin(); iHot != hots->end(); ++iHot) {
      if (iHot->isActive()) {
        double resid, residErr;
        if(iHot->resid(resid, residErr, true)){
          double value = fabs(resid/residErr);
          if (value > _maxhitchi && value > worst) {
            worst = value;
            worstHot = iHot.get();
          }
        }
      }
    }
    if(0 != worstHot){
      worstHot->setActivity(false);
      worstHot->setUsability(-5);
      retval = kalrep->fit();
      kalrep->addHistory(retval, "HitWeed");
      // Recursively iterate
      niter++;
      if (retval.success() && niter < _maxweed ) {
        retval = weedHits(kalrep, niter);
      }
    }
    return retval;
  }
}



using mu2e::KalFit;
DEFINE_FWK_MODULE(KalFit);
