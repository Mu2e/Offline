//
// Final Kalman fit
//
// Original author D. Brown and G. Tassielli
//
/*
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetector.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
// utiliites
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
// data
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/TrkReco/inc/KalFitData.hh"
#include "Offline/TrkPatRec/inc/KalFinalFit_types.hh"
#include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
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
  using namespace KalFinalFitTypes;

  class KalFinalFit : public art::EDProducer
  {
    public:
      explicit KalFinalFit(fhicl::ParameterSet const&);
      virtual ~KalFinalFit();
      void beginRun(art::Run& aRun);
    private:
      void produce(art::Event& event) override;

      unsigned _iev;
      // configuration parameters
      int _debug;
      int _diag;
      int _printfreq;
      int _cprmode;
      bool _saveall,_addhits;
      vector<double> _zsave;
      // event object tokens
      art::ProductToken<ComboHitCollection> const _shToken;
      art::InputTag const _shfTag;
      art::ProductToken<StrawHitFlagCollection> const _shfToken;
      art::ProductToken<KalSeedCollection> const _ksToken;
      art::ProductToken<KalHelixAssns> const _khaToken;
      art::ProductToken<CaloClusterCollection> const _clToken;
      // flags
      StrawHitFlag _addsel;
      StrawHitFlag _addbkg;
      TrkFitFlag _goodseed;
      double _maxdtmiss;
      // outlier cuts
      double _maxadddoca, _maxaddchi, _maxtchchi;
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      // event objects
      const ComboHitCollection* _chcol;
      const StrawHitFlagCollection* _shfcol;
      const KalSeedCollection * _kscol;
      const KalHelixAssns * _khassns;
      const CaloClusterCollection* _clCol;
      // Kalman fitter
      KalFit _kfit;
      KalFitData _result;

      // diagnostic
      Data_t                                _data;
      std::unique_ptr<ModuleHistToolBase>   _hmanager;

      // helper functions
      bool findData(const art::Event& e);
      void findMissingHits(KalFitData&kalData);
      void findMissingHits_cpr(StrawResponse::cptr_t srep, KalFitData&kalData);
      bool hasTrkCaloHit(KalFitData&kalData);

      ProditionsHandle<StrawResponse> _strawResponse_h;
      ProditionsHandle<Mu2eDetector> _mu2eDetector_h;
      ProditionsHandle<Tracker> _alignedTracker_h;
      // flow diagnostic
  };

  KalFinalFit::KalFinalFit(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _debug(pset.get<int>("debugLevel", 0)),
    _diag(pset.get<int>("diagLevel",0)),
    _printfreq(pset.get<int>("printFrequency", 101)),
    _cprmode(pset.get<int>("cprmode",0)),
    _saveall(pset.get<bool>("saveall", false)),
    _addhits(pset.get<bool>("addhits", true)),
    _zsave(pset.get<vector<double>>("ZSavePositions", vector<double>{-1522.0,0.0,1522.0})), // front, middle and back of the tracker
    _shToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _shfTag{pset.get<art::InputTag>("StrawHitFlagCollection", "none")},
    _shfToken{consumes<StrawHitFlagCollection>(_shfTag)},
    _ksToken{consumes<KalSeedCollection>(pset.get<art::InputTag>("SeedCollection"))},
    _khaToken{consumes<KalHelixAssns>(pset.get<art::InputTag>("SeedCollection"))},
    _clToken{consumes<CaloClusterCollection>(pset.get<art::InputTag>("CaloClusterCollection"))},
    _addsel(pset.get<vector<string>>("AddHitSelectionBits", vector<string>{})),
    _addbkg(pset.get<vector<string>>("AddHitBackgroundBits", vector<string>{})),
    _goodseed(pset.get<vector<string>>("GoodKalSeedFitBits", vector<string>{})),
    _maxdtmiss(pset.get<double>("DtMaxMiss",40.0)),
    _maxadddoca(pset.get<double>("MaxAddDoca",2.75)),
    _maxaddchi(pset.get<double>("MaxAddChi",4.0)),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle", TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection", TrkFitDirection::downstream))),
    _kfit(pset.get<fhicl::ParameterSet>("KalFit", {})),
    _result()
    {

      produces<KalRepCollection>();
      produces<KalRepPtrCollection>();
      //produces<StrawHitFlagCollection>();
      produces<KalSeedCollection>();
      produces<KalHelixAssns>();
      //-----------------------------------------------------------------------------
      // provide for interactive diagnostics
      //-----------------------------------------------------------------------------
      _data.result    = &_result;

      if (_diag != 0) {
        _hmanager = art::make_tool<ModuleHistToolBase>(pset.get<fhicl::ParameterSet>("diagPlugin"));
        fhicl::ParameterSet ps1 = pset.get<fhicl::ParameterSet>("KalFit.DoubletAmbigResolver");
        _data.dar               = new DoubletAmbigResolver(ps1,0,0,0);
        _data.listOfDoublets    = new std::vector<Doublet>;
        // histogram booking belongs to beginJob, KalFinalFit doesn't have it
        art::ServiceHandle<art::TFileService> tfs;
        _hmanager->bookHistograms(tfs);
      }
      else {
        _hmanager = std::make_unique<ModuleHistToolBase>();
        _data.dar            = NULL;
        _data.listOfDoublets = NULL;
      }
    }

  KalFinalFit::~KalFinalFit(){
    if (_data.dar) {
      delete _data.listOfDoublets;
      delete _data.dar;
    }
  }
  //-----------------------------------------------------------------------------
  void KalFinalFit::beginRun(art::Run& r) {
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _data.calorimeter = ch.get();

    _kfit.setCalorimeter (_data.calorimeter);
    _kfit.setCaloGeom();
  }


  void KalFinalFit::produce(art::Event& event ) {

    auto srep = _strawResponse_h.getPtr(event.id());
    auto detmodel = _mu2eDetector_h.getPtr(event.id());

    _data.tracker = _alignedTracker_h.getPtr(event.id()).get();
    _kfit.setTracker(_data.tracker);

    // event printout
    _iev=event.id().event();
    if(_debug > 0 && (_iev%_printfreq)==0)cout<<"KalFinalFit: event="<<_iev<<endl;
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::KalFinalFit: data missing or incomplete"<< endl;
    }
    // find the cluster handle (again).  This is inefficient and hard to follow FIXME!
    auto clH = event.getValidHandle(_clToken);
    // create output
    unique_ptr<KalRepCollection>    krcol(new KalRepCollection );
    unique_ptr<KalRepPtrCollection> krPtrcol(new KalRepPtrCollection );
    unique_ptr<KalSeedCollection> kscol(new KalSeedCollection());
    unique_ptr<KalHelixAssns> kfhassns (new KalHelixAssns());
    //unique_ptr<StrawHitFlagCollection> shfcol(new StrawHitFlagCollection());
    // lookup productID for payload saver
    art::ProductID kalRepsID(event.getProductID<KalRepCollection>());
    auto KalSeedCollectionPID = event.getProductID<KalSeedCollection>();
    auto KalSeedCollectionGetter = event.productGetter(KalSeedCollectionPID);
    // copy and merge hit flags
    //    size_t index(0);
    //for(auto const& ch : *_chcol) {
    //  StrawHitFlag flag(ch.flag());
    //  if(_shfcol != 0) flag.merge(_shfcol->at(index++));
    //  shfcol->push_back(flag);
    //}

    if (_diag!=0){
      _data.event  = &event;
      _data.eventNumber = event.event();
      _data.result = &_result;
      _data.tracks = krcol.get();
      _data.kscol  = kscol.get();
    }

    _result.fitType        = 1;
    _result.event          = &event ;
    _result.chcol          = _chcol ;
    //_result.shfcol         = _shfcol ;
    _result.shfcol         = nullptr; // _shfcol ;
    if (_kfit.useTrkCaloHit()) _result.caloClusterCol = _clCol;
    //    _result.tpart       = _tpart ;
    _result.fdir           = _fdir  ;

    // loop over the seed fits.  I need an index loop here to build the Ptr
    for(size_t ikseed=0; ikseed < _kscol->size(); ++ikseed) {
      KalSeed const& kseed(_kscol->at(ikseed));
      _result.kalSeed = & kseed;
      auto hptr = (*_khassns)[ikseed].second; // Ptr to the original HelixSeed
      //      _result.tpart   = kseed.particle();
      // create a Ptr for possible added CaloCluster
      art::Ptr<CaloCluster> ccPtr;
      if (kseed.caloCluster()){
        _result.caloCluster = kseed.caloCluster().get(); // should not be using KalFitData as a common block FIXME!
        ccPtr = kseed.caloCluster(); // remember the Ptr for creating the TrkCaloHitSeed and KalSeed Ptr
      }

      // only process fits which meet the requirements
      if(kseed.status().hasAllProperties(_goodseed)) {
        // check the seed has the same basic parameters as this module expects

        // if(kseed.particle() != _tpart || kseed.fitDirection() != _fdir ) {
        //   throw cet::exception("RECO")<<"mu2e::KalFinalFit: wrong particle or direction"<< endl;
        // }

        // seed should have at least 1 segment
        if(kseed.segments().size() < 1){
          throw cet::exception("RECO")<<"mu2e::KalFinalFit: no segments"<< endl;
        }
        // build a Kalman rep around this seed
        //fill the KalFitData variable
        // _result.kalSeed = &kseed;

        // _kfit.makeTrack(_shcol,kseed,krep);
        _result.init();
        _kfit.makeTrack(srep,detmodel,_result);

        // KalRep *krep = _result.stealTrack();

        if(_debug > 1){
          if(_result.krep == 0)
            cout << "No Final fit produced " << endl;
          else{
            cout << "Seed Fit HelixTraj parameters " << _result.krep->seedTrajectory()->parameters()->parameter()
              << " covariance " << _result.krep->seedTrajectory()->parameters()->covariance()
              << " NDOF = " << _result.krep->nDof()
              << " Final Fit status " << _result.krep->fitStatus()  << endl;
          }
        }
        // if successfull, try to add missing hits
        if(_addhits && _result.krep != 0 && _result.krep->fitStatus().success()){
          // first, add back the hits on this track
          //    _result.nunweediter = 0;
          _kfit.unweedHits(_result,_maxaddchi);
          if (_debug > 0) _kfit.printUtils()->printTrack(&event,_result.krep,"banner+data+hits","CalTrkFit::produce after unweedHits");

          if (_cprmode){
            findMissingHits_cpr(srep,_result);
          }else {
            findMissingHits(_result);
          }
          //check the presence of a TrkCaloHit; if it's not present, add it
          if (_kfit.useTrkCaloHit() ){
            if (!hasTrkCaloHit(_result)){
              int icc = _kfit.addTrkCaloHit(detmodel, _result);
              if(icc >=0){
                // set the CaloCluster Ptr for the TrkCaloHitSeed.
                ccPtr = art::Ptr<CaloCluster>(clH,(size_t)icc);
              }
            }
            if ( hasTrkCaloHit(_result)) _kfit.weedTrkCaloHit(_result);
            if (_diag!=0) {
              _kfit.fillTchDiag(_result);
              _data.tchDiskId  = _result.diag.diskId;
              _data.tchAdded   = _result.diag.added;
              _data.tchDepth   = _result.diag.depth;
              _data.tchDOCA    = _result.diag.doca;
              _data.tchDt      = _result.diag.dt;
              _data.tchTrkPath = _result.diag.trkPath;
              _data.tchEnergy  = _result.diag.energy;

            }
          }

          if(_result.missingHits.size() > 0){
            _kfit.addHits(srep,detmodel,_result,_maxaddchi);
          }else if (_cprmode){
            int last_iteration  = -1;
            _kfit.fitIteration(detmodel,_result,last_iteration);
          }
          if(_debug > 1)
            cout << "AddHits Fit result " << _result.krep->fitStatus()
              << " NDOF = " << _result.krep->nDof() << endl;

          //-----------------------------------------------------------------------------
          // and weed hits again to insure that addHits doesn't add junk
          //-----------------------------------------------------------------------------
          int last_iteration  = -1;
          if (_cprmode) _kfit.weedHits(_result,last_iteration);
        }
        // put successful fits into the event
        if(_result.krep != 0 && (_result.krep->fitStatus().success() || _saveall)){
          //-----------------------------------------------------------------------------
          // now evaluate the T0 and its error using the straw hits
          //-----------------------------------------------------------------------------
          //    int last_iteration  = -1;
          //    if (_cprmode) _kfit.updateT0(_result, last_iteration);

          // warning about 'fit current': this is not an error
          if(!_result.krep->fitCurrent()){
            cout << "Fit not current! " << endl;
            _result.deleteTrack();
          } else {
            // flg all hits as belonging to a track.  Doesn't work for TrkCaloHit FIXME!
            if(ikseed<StrawHitFlag::_maxTrkId){
              for(auto ihit=_result.krep->hitVector().begin();ihit != _result.krep->hitVector().end();++ihit){
                TrkStrawHit* tsh = dynamic_cast<TrkStrawHit*>(*ihit);
                //if((*ihit)->isActive() && tsh != 0)shfcol->at(tsh->index()).merge(StrawHitFlag::track);
                if (tsh == nullptr) continue;
                StrawHitFlag* flag = (StrawHitFlag*) &tsh->comboHit().flag();
                //if((*ihit)->isActive() && tsh != 0)shfcol->at(tsh->index()).merge(StrawHitFlag::track);
                if ((*ihit)->isActive()) flag->merge(StrawHitFlag::track);
              }
            }


            // save successful kalman fits in the event
            KalRep *krep = _result.stealTrack();
            krcol->push_back(krep);

            int index = krcol->size()-1;
            krPtrcol->emplace_back(kalRepsID, index, event.productGetter(kalRepsID));
            // convert successful fits into 'seeds' for persistence
            TrkFitFlag fflag(kseed.status());
            fflag.merge(TrkFitFlag::KFF);
            if(krep->fitStatus().success()) fflag.merge(TrkFitFlag::kalmanOK);
            if(krep->fitStatus().success()==1) fflag.merge(TrkFitFlag::kalmanConverged);
            //    KalSeed fseed(_tpart,_fdir,krep->t0(),krep->flt0(),kseed.status());
            KalSeed fseed(PDGCode::type(krep->particleType().particleType()),fflag,krep->flt0());
            // fill with new information
            fseed._flt0 = krep->flt0();
            // global fit information
            fseed._chisq = krep->chisq();
            fseed._ndof = krep->nDof();
            // compute the fit consistency.  Note our fit has effectively 6 parameters as t0 is allowed to float and its error is propagated to the chisquared
            fseed._fitcon =  TrkUtilities::chisqConsistency(krep);
            TrkUtilities::fillStrawHitSeeds(krep,*_chcol,fseed._hits);
            TrkUtilities::fillStraws(krep,fseed._straws);
            // sample the fit at the requested z positions.  Need options here to define a set of
            // standard points, or to sample each unique segment on the fit FIXME!
            for(auto zpos : _zsave) {
              // compute the flightlength for this z
              double fltlen = krep->pieceTraj().zFlight(zpos);
              // sample the momentum at this flight.  This belongs in a separate utility FIXME
              BbrVectorErr momerr = krep->momentumErr(fltlen);
              // sample the helix
              double locflt(0.0);
              const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(fltlen,locflt));
              // fill the segment
              KalSegment kseg;
              TrkUtilities::fillSegment(*htraj,locflt,fltlen,krep->t0(),_tpart.mass(),kseed.segments().begin()->centralHelix().charge(),_kfit.bField(),kseg);
              fseed._segments.push_back(kseg);
            }
            // see if there's a TrkCaloHit
            const TrkCaloHit* tch = TrkUtilities::findTrkCaloHit(krep);
            if(tch != 0){
              auto tmom = krep->momentum(tch->fltLen());
              TrkUtilities::fillCaloHitSeed(tch,tmom,fseed._chit);
              // set the Ptr using the helix: this could be more direct FIXME!
              fseed._chit._cluster = ccPtr;
              // create a helix segment at the TrkCaloHit
              KalSegment kseg;
              // sample the momentum at this flight.  This belongs in a separate utility FIXME
              BbrVectorErr momerr = krep->momentumErr(tch->fltLen());
              double locflt(0.0);
              const HelixTraj* htraj = dynamic_cast<const HelixTraj*>(krep->localTrajectory(tch->fltLen(),locflt));
              TrkUtilities::fillSegment(*htraj,locflt,tch->fltLen(),krep->t0(),_tpart.mass(),kseed.segments().begin()->centralHelix().charge(),_kfit.bField(),kseg);
              fseed._segments.push_back(kseg);
            }
            // save KalSeed for this track
            kscol->push_back(fseed);
            // fill assns with the helix seed
            auto kseedptr = art::Ptr<KalSeed>(KalSeedCollectionPID,kscol->size()-1,KalSeedCollectionGetter);
            kfhassns->addSingle(kseedptr,hptr);

            if (_diag > 0) _hmanager->fillHistograms(&_data);
          }
        } else {// fit failure
          _result.deleteTrack();
          //    delete krep;
        }
      }
    }

    // if (_diag > 0) _hmanager->fillHistograms(&_data);

    // put the output products into the event
    event.put(move(krcol));
    event.put(move(krPtrcol));
    event.put(move(kscol));
    event.put(move(kfhassns));
    //event.put(move(shfcol));
  }

  // find the input data objects
  bool KalFinalFit::findData(const art::Event& evt){
    _chcol = 0;
    _kscol = 0;
    _khassns = 0;

    auto shH = evt.getValidHandle(_shToken);
    _chcol = shH.product();
    auto ksH = evt.getValidHandle(_ksToken);
    _kscol = ksH.product();
    auto khaH = evt.getValidHandle(_khaToken);
    _khassns = khaH.product();
    if(_shfTag.label() != "none"){
      // auto shfH = evt.getValidHandle(_shfToken);
      //_shfcol = shfH.product();
      //if(_shfcol->size() != _chcol->size())
      //  throw cet::exception("RECO")<<"mu2e::KalFinalFit: inconsistent input collections"<< endl;
    } else {
      // _shfcol = 0;
    }
    if(_kfit.useTrkCaloHit() == 1){
      auto clH = evt.getValidHandle(_clToken);
      _clCol = clH.product();
    }

    return _chcol != 0 && _kscol != 0;
  }
  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void KalFinalFit::findMissingHits_cpr(StrawResponse::cptr_t srep, KalFitData& KRes) {

    const char* oname = "KalFinalFit::findMissingHits_cpr";

    Hep3Vector tdir;
    HepPoint   tpos;
    double     dt;

    KalRep* krep =  KRes.krep;

    krep->pieceTraj().getInfo(0.0,tpos,tdir);

    const TrkDifPieceTraj* reftraj = krep->referenceTraj();

    if (_debug > 0) printf("[%s]      shId    sec     panel       doca        drift     dr added\n",oname);

    KRes.missingHits.clear();
    //    KRes.doca.clear();

    MissingHit_t mh;

    int nstrs = KRes.chcol->size();
    for (int istr=0; istr<nstrs; ++istr) {
      mh.index = istr;
      //----------------------------------------------------------------------
      // 2015-02-11 change the selection bit for searching for missed hits
      //----------------------------------------------------------------------
      ComboHit const& sh    = _chcol->at(istr);
      if (sh.flag().hasAnyProperty(StrawHitFlag::dead)) {
        continue;
      }
      //-----------------------------------------------------------------------------
      // I think, we want to check the radial bit: if it is set, than at least one of
      // the two measured times is wrong...
      //-----------------------------------------------------------------------------
      //      int radius_ok = _shfcol->at(istr).hasAllProperties(StrawHitFlag::radsel);
      dt        = _chcol->at(istr).time()-KRes.krep->t0()._t0;

      //      if (radius_ok && (fabs(dt) < _maxdtmiss)) {
      if (fabs(dt) < _maxdtmiss) {
        // make sure we haven't already used this hit

        TrkStrawHit  *tsh, *closest(NULL);
        bool found = false;

        Straw const&      straw = _data.tracker->getStraw(sh.strawId());
        CLHEP::Hep3Vector hpos  = straw.getMidPoint();

        double            dz_max(1.e12) ; // closest_z(1.e12);
        double            zhit = hpos.z();

        for (std::vector<TrkHit*>::iterator it=krep->hitVector().begin(); it!=krep->hitVector().end(); it++) {
          //    tsh = static_cast<TrkStrawHit*> (*it);
          tsh = dynamic_cast<TrkStrawHit*> (*it);
          if (tsh ==0)                  continue;
          int tsh_index = tsh->index();
          if (tsh_index == istr) {
            found = true;
            break;
          }
          // check proximity in Z
          Straw const&  trk_straw = _data.tracker->getStraw(tsh->comboHit().strawId());
          double        ztrk      = trk_straw.getMidPoint().z();

          double dz  = ztrk-zhit;
          if (fabs(dz) < fabs(dz_max)) {
            closest   = tsh;
            dz_max    = dz;
          }
        }

        if (! found) {
          // estimate trajectory length to hit
          double hflt = 0;
          TrkHelixUtils::findZFltlen(*reftraj,zhit,hflt);

          // good in-time hit.  Compute DOCA of the wire to the trajectory
          // also estimate the drift time if this hit were on the track to get hte hit radius

          TrkT0 hitt0 = closest->hitT0();

          double s    = closest->fltLen();
          double mom  = krep->momentum(s).mag();
          double beta = krep->particleType().beta(mom);

          double tflt = (hflt-s)/(beta*CLHEP::c_light);
          hitt0._t0  += tflt;

          CLHEP::Hep3Vector hdir  = straw.getDirection();
          // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
          HepPoint          spt(hpos.x(),hpos.y(),hpos.z());
          TrkLineTraj       htraj(spt,hdir,-20,20);
          // estimate flightlength along track.  This assumes a constant BField!!!

          double      fltlen = (hpos.z()-tpos.z())/tdir.z();

          TrkPoca     hitpoca(krep->pieceTraj(),fltlen,htraj,0.0);

          double      rdrift;//, hit_error(0.2);

          TrkStrawHit hit(srep,sh,*_data.tracker,istr,hitt0,hflt,1.,1.);//hit_error,1.,_maxadddoca,1.);

          double tdrift=hit.time()-hit.hitT0()._t0;

          //    tcal->TimeToDistance(straw.index(),tdrift,tdir,t2d);

          // find the track direction at this hit
          Hep3Vector tdir  = krep->traj().direction(fltlen);
          Hep3Vector tperp = tdir - tdir.dot(straw.getDirection())*straw.getDirection();
          double     phi   = tperp.theta();

          rdrift = srep->driftTimeToDistance(straw.id(),tdrift,phi);

          mh.doca   = hitpoca.doca();
          if (mh.doca > 0) mh.dr = mh.doca-rdrift;
          else             mh.dr = mh.doca+rdrift;
          //-----------------------------------------------------------------------------
          // flag hits with small residuals
          //-----------------------------------------------------------------------------
          int added (0);
          if (fabs(mh.dr) < _maxadddoca) {
            KRes.missingHits.push_back(mh);
            //            KRes.doca.push_back(doca);
            added = 1;
          }

          if (_debug > 0) {
            printf("[CalTrkFit::findMissingHits] %8i  %6i  %8i  %10.3f %10.3f %10.3f   %3i\n",
                straw.id().asUint16(),
                straw.id().getPlane(),
                straw.id().getPanel(),
                mh.doca, rdrift, mh.dr,added);
          }
        }
      }
      else {
        if (_debug > 0) {
          printf("[%s] rejected hit: i, index, flag, dt: %5i %5i %s %10.3f\n",
              oname,istr,sh.strawId().asUint16(),
              KRes.shfcol->at(istr).hex().data(),-1.);
        }
      }
    }
    //-----------------------------------------------------------------------------
    // sort hits in accending |dr|, such that the first hit has the smallest |dr|
    //-----------------------------------------------------------------------------
    int nmiss = KRes.missingHits.size();
    for (int i1=0; i1<nmiss-1; i1++) {
      MissingHit_t* h1 = &KRes.missingHits[i1];
      for (int i2=i1+1; i2<nmiss; i2++) {
        MissingHit_t* h2 = &KRes.missingHits[i2];
        if (fabs(h1->dr) > fabs(h2->dr)) {
          // swap hits
          MissingHit_t tmp = *h1;
          *h1 = *h2;
          *h2 = tmp;
        }
      }
    }
    }

    void KalFinalFit::findMissingHits(KalFitData&kalData) {
      KalRep* krep = kalData.krep;

      //clear the array
      kalData.missingHits.clear();
      const Tracker& tracker = *_data.tracker;

      //  Trajectory info
      Hep3Vector tdir;
      HepPoint tpos;
      krep->pieceTraj().getInfo(krep->flt0(),tpos,tdir);
      unsigned nstrs = _chcol->size();
      TrkStrawHitVector tshv;
      convert(krep->hitVector(),tshv);
      for(unsigned istr=0; istr<nstrs;++istr){
        //if(_shfcol->at(istr).hasAllProperties(_addsel)&& !_shfcol->at(istr).hasAnyProperty(_addbkg)){
        //  ComboHit const& sh = _chcol->at(istr);
        //  if (sh.flag().hasAnyProperty(StrawHitFlag::dead)) {
         ComboHit const& sh = _chcol->at(istr);
         StrawHitFlag flag = sh.flag();

         // if(_shfcol->at(istr).hasAllProperties(_addsel)&& !_shfcol->at(istr).hasAnyProperty(_addbkg)){
         if(flag.hasAllProperties(_addsel)&& !flag.hasAnyProperty(_addbkg)){
           if (flag.hasAnyProperty(StrawHitFlag::dead)) {
             continue;
          }
          if(fabs(_chcol->at(istr).time()-krep->t0()._t0) < _maxdtmiss) {
            // make sure we haven't already used this hit
            vector<TrkStrawHit*>::iterator ifnd = find_if(tshv.begin(),tshv.end(),FindTrkStrawHit(sh));
            if(ifnd == tshv.end()){
              // good in-time hit.  Compute DOCA of the wire to the trajectory
              Straw const& straw = tracker.getStraw(sh.strawId());
              CLHEP::Hep3Vector hpos = straw.getMidPoint();
              CLHEP::Hep3Vector hdir = straw.getDirection();
              // convert to HepPoint to satisfy antique BaBar interface: FIXME!!!
              HepPoint spt(hpos.x(),hpos.y(),hpos.z());
              TrkLineTraj htraj(spt,hdir,-straw.halfLength(),straw.halfLength());
              // estimate flightlength along track.  This assumes a constant BField!!!
              double fltlen = (hpos.z()-tpos.z())/tdir.z();
              // estimate hit length
              HepPoint tp = krep->pieceTraj().position(fltlen);
              Hep3Vector tpos(tp.x(),tp.y(),tp.z()); // ugly conversion FIXME!
              double hitlen = hdir.dot(tpos - hpos);
              TrkPoca hitpoca(krep->pieceTraj(),fltlen,htraj,hitlen);

              // flag hits with small residuals
              if(fabs(hitpoca.doca()) < _maxadddoca){
                MissingHit_t m;
                m.index = istr;
                m.doca  = hitpoca.doca();
                // m.dr = ??;
                kalData.missingHits.push_back(m);
              }
            }
          }
        }
      }
    }

    //--------------------------------------------------------------------------------
    // function to check the presence of a TrkCaloHit in the KalRep
    //--------------------------------------------------------------------------------
    bool KalFinalFit::hasTrkCaloHit(KalFitData&kalData){
      bool retval(false);

      TrkHitVector *thv      = &(kalData.krep->hitVector());
      for (auto ihit=thv->begin();ihit!=thv->end(); ++ihit){
        TrkCaloHit*hit = dynamic_cast<TrkCaloHit*>(*ihit);
        if (hit != 0){
          retval = true;
          break;
        }
      }

      return retval;
    }


  }// mu2e

  DEFINE_ART_MODULE(mu2e::KalFinalFit)
*/
