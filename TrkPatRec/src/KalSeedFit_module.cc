//
// Starting with a Helix fit, perform a Least-squares fit (using the BTrk
// Kalman fit, appropriately configured) to produce an initial estimate of the parmeters and
// covariance for the final Kalman fit.  This fit uses wire positions only,
// not drift.
//
// Original author Dave Brown (LBNL) 31 Aug 2016
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Run.h"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/Mu2eMaterial.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetector.hh"
// utiliites
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/TrkReco/inc/TrkDef.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
// data
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
// Diagnostic data objects
#include "Offline/TrkReco/inc/KalFitData.hh"
#include "Offline/TrkPatRec/inc/KalSeedFit_types.hh"
// BaBar
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
// Mu2e BaBar
#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "Offline/TrkReco/inc/KalFit.hh"
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

  using namespace KalSeedFitTypes;


  class KalSeedFit : public art::EDProducer
  {
    public:
      explicit KalSeedFit(fhicl::ParameterSet const&);
      virtual ~KalSeedFit();
      virtual void beginRun(art::Run&);
      virtual void produce(art::Event& event );
    private:
      unsigned _iev;
      // configuration parameters
      int _debug;
      int _diag;
      int _printfreq;
      bool _saveall;
      bool _checkhelicity;
      // event object tags
      art::ProductToken<ComboHitCollection> const _shToken;
      art::ProductToken<HelixSeedCollection> const _hsToken;
      TrkFitFlag _seedflag; // helix fit flag
      unsigned _minnhits; // minimum # of hits
      double _maxdoca;      // outlier cut
      bool _foutliers; // filter hits far from the helix
      bool _fhoutliers; // filter hits found flagged as outliers in the helix fit
      double _maxAddDoca;   // rescue hits cut after fit
      double _maxAddChi;    // cut for KalFit::AddHits
      int    _rescueHits;   // search for missing hits after the fit is performed
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      vector<double> _perr; // diagonal parameter errors to use in the fit
      Helicity _helicity; // cached value of helicity expected for this fit
      double _upz, _downz; // z positions to extend the segment
      double _amsign; // cached sign of angular momentum WRT the z axis
      double _bz000;        // sign of the magnetic field at (0,0,0)
      HepSymMatrix _hcovar; // cache of parameter error covariance matrix
      TrkFitFlag  _ksf; // default fit flag
      // cache of event objects
      const ComboHitCollection *_chcol;
      const HelixSeedCollection *_hscol;
      // ouptut collections
      // Kalman fitter.  This will be configured for a least-squares fit (no material or BField corrections).
      KalFit _kfit;
      KalFitData _result;
      const Tracker* _tracker;     // straw tracker geometry

      ProditionsHandle<StrawResponse> _strawResponse_h;
      ProditionsHandle<Mu2eMaterial> _mu2eMaterial_h;
      ProditionsHandle<Mu2eDetector> _mu2eDetector_h;
      ProditionsHandle<Tracker> _alignedTracker_h;

      // diagnostic
      Data_t                                _data;
      std::unique_ptr<ModuleHistToolBase>   _hmanager;

      // helper functions
      bool findData(const art::Event& e);
      void filterOutliers(TrkDef& trkdef);
      void findMissingHits(KalFitData&kalData);
  };

  KalSeedFit::KalSeedFit(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _debug(pset.get<int>("debugLevel",0)),
    _diag(pset.get<int>("diagLevel",0)),
    _printfreq(pset.get<int>("printFrequency",101)),
    _saveall(pset.get<bool>("saveall",false)),
    _checkhelicity(pset.get<bool>("CheckHelicity",true)),
    _shToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _hsToken{consumes<HelixSeedCollection>(pset.get<art::InputTag>("SeedCollection"))},
    _seedflag(pset.get<vector<string> >("HelixFitFlag",vector<string>{"HelixOK"})),
    _minnhits(pset.get<unsigned>("MinNHits",10)),
    _maxdoca(pset.get<double>("MaxDoca",40.0)),
    _foutliers(pset.get<bool>("FilterOutliers",true)),
    _fhoutliers(pset.get<bool>("FilterHelixOutliers",false)),
    _maxAddDoca(pset.get<double>("MaxAddDoca")),
    _maxAddChi(pset.get<double>("MaxAddChi")),
    _rescueHits(pset.get<int>("rescueHits")),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _perr(pset.get<vector<double> >("ParameterErrors")),
    _upz(pset.get<double>("UpstreamZ",-1500)),
    _downz(pset.get<double>("DownstreamZ",1500)),
    _ksf(TrkFitFlag::KSF),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit",fhicl::ParameterSet())),
    _result()
    {
      // This following consumesMany call is necessary because
      // ComboHitCollection::fillStrawHitIndices calls getManyByType
      // under the covers.
      consumesMany<ComboHitCollection>();
      produces<KalSeedCollection>();
      produces<KalHelixAssns>();
      // check dimensions
      if(_perr.size() != HelixTraj::NHLXPRM)
        throw cet::exception("RECO")<<"mu2e::KalSeedFit: parameter error vector has wrong size"<< endl;
      // mock covariance matrix, all diagonal
      _hcovar = HepSymMatrix(HelixTraj::NHLXPRM,0);
      for(size_t ipar = 0; ipar < HelixTraj::NHLXPRM; ++ipar){
        _hcovar(ipar+1,ipar+1) = _perr[ipar]*_perr[ipar]; // clhep indexing starts a 1
      }

      //-----------------------------------------------------------------------------
      // provide for interactive disanostics
      //-----------------------------------------------------------------------------
      _data.result    = &_result;


      if (_diag != 0) _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
      else            _hmanager = std::make_unique<ModuleHistToolBase>();
    }

  KalSeedFit::~KalSeedFit(){}

  void KalSeedFit::beginRun(art::Run& run){
    // calculate the helicity
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;

    // initialize the BTrk material and particle models
    _mu2eMaterial_h.get(run.id());

    _kfit.setCaloGeom();

    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(0.0,0.0,0.0);
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    // helicity is a purely geometric quantity, however it's easiest
    // to determine it from the kinematics (angular momentum and Z momentum)
    _bz000    = field.z();
    _amsign   = copysign(1.0,-_tpart.charge()*_bz000);
    _helicity = Helicity(static_cast<float>(_fdir.dzdt()*_amsign));
  }

  void KalSeedFit::produce(art::Event& event ) {

    auto srep = _strawResponse_h.getPtr(event.id());
    auto detmodel = _mu2eDetector_h.getPtr(event.id());

    _tracker = _alignedTracker_h.getPtr(event.id()).get();
    _kfit.setTracker(_tracker);

    // create output collection
    unique_ptr<KalSeedCollection> kscol(new KalSeedCollection());
    unique_ptr<KalHelixAssns> ksha(new KalHelixAssns());
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"KalSeedFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::KalSeedFit: data missing or incomplete"<< endl;
    }
    if (_diag){
      _data.event  = &event;
      _data.result = &_result;
      _data.nrescued.clear();
      // _data.mom.clear();
      _data.tracks = kscol.get();
    }

    _result.fitType     = 0;
    _result.event       = &event ;
    _result.chcol       = _chcol ;
    //    _result.tpart       = _tpart ;
    _result.fdir        = _fdir  ;

    // loop over the Helices
    for (size_t iseed=0; iseed<_hscol->size(); ++iseed) {
      // convert the HelixSeed to a TrkDef
      HelixSeed const& hseed(_hscol->at(iseed));

      if (hseed.caloCluster()) _result.caloCluster = hseed.caloCluster().get();
      _result.helixSeed = &hseed;
      //-----------------------------------------------------------------------------
      // 2018-12-08 PM : allow list of helices to contain helices of different
      // helicities and corresponding to particles of opposite signs. Assume that the
      // PDG particle coding scheme is used such that the particle and antiparticle
      // PDG codes have opposite signs
      //-----------------------------------------------------------------------------
      TrkParticle tpart(_tpart);
      if(_helicity != hseed.helix().helicity()) {
        if(_checkhelicity) throw cet::exception("RECO")<<"mu2e::KalSeedFit: helicity doesn't match configuration" << endl;
        TrkParticle::type t = (TrkParticle::type) (-int(_tpart.particleType()));
        tpart = TrkParticle(t);
      }

      double amsign   = copysign(1.0,-tpart.charge()*_bz000);

      HepVector hpvec(HelixTraj::NHLXPRM);
      // verify the fit meets requirements and can be translated
      // to a fit trajectory.  This accounts for the physical particle direction
      // helicity.  This could be wrong due to FP effects, so don't treat it as an exception
      if(hseed.status().hasAllProperties(_seedflag) &&
          //   _helicity == hseed.helix().helicity() &&
          TrkUtilities::RobustHelix2Traj(hseed._helix,hpvec,amsign)){
        HelixTraj hstraj(hpvec,_hcovar);
        // update the covariance matrix
        if(_debug > 1)
          //    hstraj.printAll(cout);
          cout << "Seed Fit HelixTraj parameters " << hstraj.parameters()->parameter()
            << "and covariance " << hstraj.parameters()->covariance() <<  endl;
        // build a time cluster: exclude the outlier hits
        TimeCluster tclust;
        tclust._t0 = hseed._t0;
        for(uint16_t ihit=0;ihit < hseed.hits().size(); ++ihit){
          ComboHit const& ch = hseed.hits()[ihit];
          if((!_fhoutliers) || (!ch.flag().hasAnyProperty(StrawHitFlag::outlier)))
            hseed.hits().fillStrawHitIndices(ihit,tclust._strawHitIdxs);
        }
        // create a TrkDef; it should be possible to build a fit from the helix seed directly FIXME!
        //  TrkDef seeddef(tclust,hstraj,_tpart,_fdir);
        TrkDef seeddef(tclust,hstraj,tpart,_fdir);
        // filter outliers; this doesn't use drift information, just straw positions
        if(_foutliers)filterOutliers(seeddef);
        const HelixTraj* htraj = &seeddef.helix();
        double           flt0  = htraj->zFlight(0.0);
        double           mom   = TrkMomCalculator::vecMom(*htraj, _kfit.bField(), flt0).mag();
        double           vflt  = seeddef.particle().beta(mom)*CLHEP::c_light;
        double           helt0 = hseed.t0().t0();

        KalSeed kf(PDGCode::type(tpart.particleType()),hseed.status(), flt0 );
        // extract the hits from the rep and put the hitseeds into the KalSeed
        int nsh = seeddef.strawHitIndices().size();//tclust._strawHitIdxs.size();
        for (int i=0; i< nsh; ++i){
          size_t          istraw   = seeddef.strawHitIndices().at(i);
          const ComboHit& strawhit(_chcol->at(istraw));
          const Straw&    straw    = _tracker->getStraw(strawhit.strawId());
          double          fltlen   = htraj->zFlight(straw.getMidPoint().z());
          double          propTime = (fltlen-flt0)/vflt;

          //fill the TrkStrwaHitSeed info
          TrkStrawHitSeed tshs;
          tshs._index  = istraw;
          tshs._t0     = TrkT0(helt0 + propTime, hseed.t0().t0Err());
          tshs._trklen = fltlen;
          kf._hits.push_back(tshs);
        }

        if(kf._hits.size() >= _minnhits) kf._status.merge(TrkFitFlag::hitsOK);
        // extract the helix trajectory from the fit (there is just 1)
        // use this to create segment.  This will be the only segment in this track
        if(htraj != 0){
          KalSegment kseg;
          // sample the momentum at this point
          TrkUtilities::fillSegment(*htraj,0.0,0.0,hseed.t0(),tpart.mass(),int(tpart.charge()),_kfit.bField(),kseg);
          kf._segments.push_back(kseg);
        } else {
          throw cet::exception("RECO")<<"mu2e::KalSeedFit: Can't extract helix traj from seed fit" << endl;
        }

        // now, fit the seed helix from the filtered hits

        //fill the KalFitData variable
        _result.kalSeed = &kf;

        _kfit.makeTrack(srep,detmodel,_result);

        if(_debug > 1){
          if(_result.krep == 0)
            cout << "No Seed fit produced " << endl;
          else
            cout << "Seed Fit result " << _result.krep->fitStatus()  << endl;
        }
        if(_result.krep != 0 && (_result.krep->fitStatus().success() || _saveall)){
          if (_rescueHits) {
            int nrescued = 0;
            findMissingHits(_result);
            nrescued = _result.missingHits.size();
            if (nrescued > 0) {
              _kfit.addHits(srep,detmodel,_result, _maxAddChi);
            }
          }

          //    KalRep *krep = _result.stealTrack();

          // convert the status into a FitFlag
          // create a KalSeed object from this fit, recording the particle and fit direction
          //    KalSeed kseed(_tpart,_fdir,_result.krep->t0(),_result.krep->flt0(),seedok);

          KalSeed kseed(PDGCode::type(_result.krep->particleType().particleType()),kf.status(), _result.krep->flt0());
          kseed._status.merge(_ksf);
          if(_result.krep->fitStatus().success())kseed._status.merge(TrkFitFlag::kalmanOK);
          // add CaloCluster if present
          kseed._chit._cluster = hseed.caloCluster();
          // extract the hits from the rep and put the hitseeds into the KalSeed
          TrkUtilities::fillStrawHitSeeds(_result.krep,*_chcol,kseed._hits);
          if(_result.krep->fitStatus().success())kseed._status.merge(TrkFitFlag::seedOK);
          if(_result.krep->fitStatus().success()==1)kseed._status.merge(TrkFitFlag::seedConverged);
          if(kseed._hits.size() >= _minnhits)kseed._status.merge(TrkFitFlag::hitsOK);
          kseed._chisq = _result.krep->chisq();
          kseed._ndof = _result.krep->nDof();
          // use the default consistency calculation, as t0 is not fit here
          kseed._fitcon = _result.krep->chisqConsistency().significanceLevel();
          // extract the helix trajectory from the fit (there is just 1)
          double locflt;
          const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(_result.krep->localTrajectory(_result.krep->flt0(),locflt));
          // use this to create segment.  This will be the only segment in this track
          if(htraj != 0){
            KalSegment kseg;
            // sample the momentum at this point
            BbrVectorErr momerr = _result.krep->momentumErr(_result.krep->flt0());
            TrkUtilities::fillSegment(*htraj,locflt,_result.krep->flt0(),_result.krep->t0(),tpart.mass(),int(tpart.charge()),_kfit.bField(),kseg);
            // extend the segment
            double upflt(0.0), downflt(0.0);
            TrkHelixUtils::findZFltlen(*htraj,_upz,upflt);
            TrkHelixUtils::findZFltlen(*htraj,_downz,downflt);
            double tup = kseg.fltToTime(upflt);
            double tdown = kseg.fltToTime(downflt);
            if(_fdir == TrkFitDirection::downstream){
              kseg._tmin = tup;
              kseg._tmax = tdown;
            } else {
              kseg._tmax = tup;
              kseg._tmin = tdown;
            }
            kseed._segments.push_back(kseg);
            // push this seed into the collection
            kscol->push_back(kseed);
            // fill assns with the helix seed
            auto hsH = event.getValidHandle(_hsToken);
            auto hptr = art::Ptr<HelixSeed>(hsH,iseed);
            auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kscol->size()-1,KalSeedCollectionGetter);
            ksha->addSingle(kseedptr,hptr);
            if(_debug > 1){
              cout << "Seed fit segment parameters " << endl;
              for(size_t ipar=0;ipar<5;++ipar) cout << kseg.helix()._pars[ipar] << " ";
              cout << " covariance " << endl;
              for(size_t ipar=0;ipar<15;++ipar)
                cout << kseg.covar()._cov[ipar] << " ";
              cout << endl;
            }
          } else {
            throw cet::exception("RECO")<<"mu2e::KalSeedFit: Can't extract helix traj from seed fit" << endl;
          }
        }
        // cleanup the seed fit KalRep.  Optimally the krep should be a data member of this module
        // and get reused to avoid thrashing memory, but the BTrk code doesn't support that, FIXME!
        _result.deleteTrack();
      }
    }
    // put the tracks into the event
    event.put(move(kscol));
    event.put(move(ksha));
  }


  // find the input data objects
  bool KalSeedFit::findData(const art::Event& evt){
    _chcol = 0;
    _hscol = 0;

    auto shH = evt.getValidHandle(_shToken);
    _chcol = shH.product();
    auto hsH = evt.getValidHandle(_hsToken);
    _hscol = hsH.product();

    return _chcol != 0 && _hscol != 0;
  }

  void KalSeedFit::filterOutliers(TrkDef& mydef){
    // for now filter on DOCA.  In future this shoudl be an MVA using time and position FIXME!
    //  Trajectory info
    Hep3Vector tdir;
    HepPoint tposp;
    double flt0 = mydef.helix().zFlight(0.0);
    mydef.helix().getInfo(flt0,tposp,tdir);
    // tracker and conditions
    const Tracker& tracker = *_tracker;


    const vector<StrawHitIndex>& indices = mydef.strawHitIndices();
    vector<StrawHitIndex> goodhits;
    for(unsigned ihit=0;ihit<indices.size();++ihit){
      ComboHit const& sh = _chcol->at(indices[ihit]);
      Straw const& straw = tracker.getStraw(sh.strawId());
      CLHEP::Hep3Vector hpos = straw.getMidPoint();
      CLHEP::Hep3Vector hdir = straw.getDirection();
      // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
      HepPoint spt(hpos.x(),hpos.y(),hpos.z());
      TrkLineTraj htraj(spt,hdir,-straw.halfLength(),straw.halfLength());
      // estimate flightlength along track.  This assumes a constant BField!!!
      double fltlen = (hpos.z()-tposp.z())/tdir.z();
      HepPoint tp = mydef.helix().position(fltlen);
      Hep3Vector tpos(tp.x(),tp.y(),tp.z()); // ugly conversion FIXME!
      double hitlen = hdir.dot(tpos - hpos);
      TrkPoca hitpoca(mydef.helix(),fltlen,htraj,hitlen);
      // keep hits with small residuals
      if(fabs(hitpoca.doca()) < _maxdoca){
        goodhits.push_back(indices[ihit]);
      }
    }
    // update track
    mydef.strawHitIndices() = goodhits;
  }

  //-----------------------------------------------------------------------------
  // look for hits which were not a part of the helix hit list around the
  // trajectory found by the seed fit
  // look at all hits included into the corresponding time cluster
  // first reactivate already associated hits
  //-----------------------------------------------------------------------------
  void KalSeedFit::findMissingHits(KalFitData&kalData){

    const char* oname = "KalSeedFit::findMissingHits";

    mu2e::TrkStrawHit*       hit;
    int                      hit_index;
    const ComboHit*          sh;
    const Straw*             straw;

    Hep3Vector               tdir;
    HepPoint                 tpos;
    double                   doca, /*rdrift, */fltlen;

    if (_debug > 0) printf("[%s]: BEGIN\n",oname);

    const KalRep* krep = kalData.krep;

    kalData.missingHits.clear();

    const TrkDifTraj& trajectory = krep->traj();
    const vector<TrkHit*>&  trackHits  = krep->hitVector();
    //-----------------------------------------------------------------------------
    // get track position and direction at S=0
    //-----------------------------------------------------------------------------
    trajectory.getInfo(0.0,tpos,tdir);
    //-----------------------------------------------------------------------------
    // look for so far unused hits around the trajectory
    //-----------------------------------------------------------------------------
    const HelixSeed*   hseed = kalData.helixSeed;
    const  std::vector<StrawHitIndex>& tchits = hseed->timeCluster()->hits();

    int n = tchits.size();
    for (int i=0; i<n; ++i) {
      hit_index = tchits.at(i);
      sh        = &kalData.chcol->at(hit_index);
      if (sh->flag().hasAnyProperty(StrawHitFlag::dead)) {
        continue;
      }
      straw     = &_tracker->getStraw(sh->strawId());

      const CLHEP::Hep3Vector& wpos = straw->getMidPoint();
      const CLHEP::Hep3Vector& wdir = straw->getDirection();

      HepPoint      wpt  (wpos.x(),wpos.y(),wpos.z());
      TrkLineTraj   wire (wpt,wdir,-20,20);
      //-----------------------------------------------------------------------------
      // estimate flightlength along the track for z-coordinate corresponding to the
      // wire position. This assumes a constant BField!!!
      // in principle, this should work well enough, however, may want to check
      // then determine the distance from the wire to the trajectory
      //-----------------------------------------------------------------------------
      fltlen = (wpos.z()-tpos.z())/tdir.z();
      TrkPoca wpoca(trajectory,fltlen,wire,0.0);
      doca   = wpoca.doca();

      int found(-1);

      if (std::fabs(doca) < _maxAddDoca) {
        found = 0;
        for (auto it=trackHits.begin(); it<trackHits.end(); it++) {
          hit    = dynamic_cast<mu2e::TrkStrawHit*> (*it);
          if (hit == 0)                   continue;     //it means that "hit" is a TrkCaloHit
          int shIndex = int(hit->index());
          if (hit_index == shIndex) {
            found = 1;
            break;
          }
        }
        //-----------------------------------------------------------------------------
        // KalSeedFit doesn't look at the hit residuals, only wires
        //-----------------------------------------------------------------------------
        if (found == 0) {
          MissingHit_t m;
          m.index = hit_index;
          m.doca  = doca;
          //    m.dr    = ??;
          kalData.missingHits.push_back(m);//hit_index);
          //    KRes._doca.push_back(doca);
        }
      }

      if (_debug > 0) printf("[%s] %5i %8.3f %2i \n",oname,hit_index,doca,found);

    }
  }

}// mu2e
using mu2e::KalSeedFit;
DEFINE_ART_MODULE(KalSeedFit)
