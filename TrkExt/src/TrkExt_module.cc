//
//  Track Extrapolation module.
//  All in Detector (=Tracker) coordinate. Points other than tracker coordinates are properly transformed. 
//  B field basically requires mu2e coordinate
//  G4 uses G4World coordinate
//
//  $Id: TrkExt_module.cc,v 1.2 2012/08/18 22:27:56 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/18 22:27:56 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/WorldG4.hh"

using namespace CLHEP;
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TH1F.h"
#include "TTree.h"


#include "KalmanTests/inc/KalRepCollection.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/PointTrajectory.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "TrkBase/TrkParticle.hh"

#include "BFieldGeom/inc/BFieldManager.hh"

#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "TrkExt/inc/TrkExtDetectors.hh"

using namespace std;

namespace mu2e {

  const double MEC22 = 0.510998910 * 0.510998910;
  const double VELOCITY_OF_LIGHT = 2.99792458e8; 
  const int MAXSIM = 5000;
  const int MAXNBACK = 5000;
  const double RUNGE_KUTTA_KQ = -1.e-9*VELOCITY_OF_LIGHT; //k = 2.99e-1, q = -1


  namespace TrkExtExitCode {
    enum Enum {
      Undefined = -1,
      DoNotExit = 0,
      WriteData = 1,
      UndefinedVolume = 2,
      MaximumMomentum = 3,
      ReflectionLimit = 4,
      MaximumPoints = 5
    };
  };


  class TrkExt : public art::EDProducer {

  public:
    explicit TrkExt(fhicl::ParameterSet const& pset);
    virtual ~TrkExt() { }
    void beginJob();
    void beginRun(art::Run &run);
    void beginSubRun(art::SubRun & lblock );
    virtual void produce(art::Event& event);
    void endJob();

  private:

    std::string _fitterModuleLabel;
    std::string _g4ModuleLabel;
    std::string _makerModuleLabel;
    double _maxMomentum;
    bool _turnOnMaterialEffect;
    bool _useStoppingPower;
    int _maxNBack;
    double _extrapolationStep; //in mm
    double _recordingStep;
    bool _mcFlag;
    bool _useVirtualDetector;
    int _bFieldGradientMode;
    bool _turnOnMultipleScattering;
    int _debugLevel;

    bool _flagEloss;
    bool _flagDiagnostics;

    TTree * _hEloss;
    TH1F * _hExitCode;
    TH1F * _hNSteps;
    TH1F * _hNData;
    TH1F * _hFL;
    TH1F * _hNPAClust;
    TH1F * _hNSTClust;
    TH1F * _hPinit;
    TH1F * _hPfin;

    float _vdx[5], _vdy[5], _vdz[5], _vdpx[5], _vdpy[5], _vdpz[5], _vdp[5];
    int _evtid, _trkid;
    float _eloss_mean, _eloss_mp, _eloss_p;
    int _eloss_vid;

    TrkExtTraj _traj;

    Hep3Vector _origin;
    Hep3Vector _mu2eOriginInWorld;

    BFieldManager const * _bfMgr;
    TrkExtDetectors _mydet;
    string _trkPatRecInstanceName;


    double _dummyStoppingTarget_halfLength;
    double _dummyStoppingTarget_z0;

    void readTrkPatRec(KalRep const & trk, 
                      Hep3Vector * xstart, Hep3Vector * pstart, 
                      Hep3Vector * xstop, Hep3Vector * pstop, 
                      HepMatrix * covstart, HepMatrix * covstop) ;
    bool readVD (const art::Event& event, TrkHotList const * hits) ;
    int doExtrapolation (Hep3Vector x, Hep3Vector p, HepMatrix cov, bool direction) ;

    HepVector _runge_kutta_newpar_5th (HepVector r0, double ds, bool mode) ;
    HepVector _runge_kutta_newpar_f (HepVector r, Hep3Vector B) ;

    TrkExtTrajPoint calculateNextPosition(TrkExtTrajPoint r00, double ds);

    Hep3Vector getBField (Hep3Vector& x) ; // in Detector coordinate
    Hep3Vector getBField (HepVector& r) ; // in Detector coordinate
    Hep3Vector getBFieldWithGradient( Hep3Vector & x, 
                                      double & bxx, double & bxy, double & bxz, 
                                      double & byx, double & byy, double & byz, 
                                      double & bzx, double & bzy, double & bzz);
    bool checkOutofReflectionLimit (Hep3Vector & x, Hep3Vector & p); // in Detector coordinate
    HepMatrix getCovarianceTransport(TrkExtTrajPoint & r0, double ds, double deltapp);
    HepMatrix getCovarianceMultipleScattering(TrkExtTrajPoint & r0, double ds);


  };

  TrkExt::TrkExt(fhicl::ParameterSet const& pset):
    _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel")),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel")),
    _maxMomentum(pset.get<double>("maxMomentum", 104.96)),
    _turnOnMaterialEffect(pset.get<bool>("turnOnMaterialEffect", true)),
    _useStoppingPower(pset.get<bool>("useStoppingPower", false)),
    _maxNBack(pset.get<int>("maxNBack", 5000)),
    _extrapolationStep(pset.get<double>("extrapolationStep", 5.0)),    // in mm
    _recordingStep(pset.get<double>("recordingStep", 10.0)),    // in mm
    _mcFlag(pset.get<bool>("mcFlag", false)),
    _useVirtualDetector(pset.get<bool>("useVirtualDetector", false)),
    _bFieldGradientMode(pset.get<int>("bFieldGradientMode", 1)),
    _turnOnMultipleScattering(pset.get<bool>("turnOnMultipleScattering", true)),
    _debugLevel(pset.get<int>("debugLevel", 1)),
    _hEloss(0)
  {

    if (_g4ModuleLabel != "") {
      _mcFlag = true;
    }

    _dummyStoppingTarget_halfLength = pset.get<double>("dummyStoppingTarget.halfLength", 400.);
    _dummyStoppingTarget_z0 = pset.get<double>("dummyStoppingTarget.z0", 5900.);

    _trkPatRecInstanceName = TrkFitDirection(TrkFitDirection::downstream).name() + TrkParticle(TrkParticle::e_minus).name();

    switch (_debugLevel) {
      case 1:
        _flagEloss = false;
        _flagDiagnostics = true;
        break;
      case 2:
        _flagEloss = true;
        _flagDiagnostics = true;
        break;
      case 0:
      default:
        _flagEloss = false;
        _flagDiagnostics = false;
        break;
    }

    produces<TrkExtTrajCollection>(); 

  }

  void TrkExt::beginJob(){

    // some input check
    if (_maxNBack > MAXNBACK) _maxNBack = MAXNBACK;
    if (_extrapolationStep == 0) _extrapolationStep = 5.;
    else if (_extrapolationStep <0) _extrapolationStep = fabs(_extrapolationStep);
    if (_recordingStep <0) _recordingStep = 0.0;
    if (!_mcFlag) {
      if (_useVirtualDetector) {
        cerr << "TrkExt: VirtualDetector turned off for data" << endl;
        _useVirtualDetector = false;
      }
    }

    if (_bFieldGradientMode != 0 
        && _bFieldGradientMode != 1) {
      cerr << "TrkExt: bFieldGradientMode forced to 1" << endl;
      _bFieldGradientMode = 1;
    }

    cerr << "TrkExt: extrapolationStep = " << _extrapolationStep << endl;
    cerr << "TrkExt: recordingStep = " << _recordingStep << endl;

    // histograms

    art::ServiceHandle<art::TFileService> tfs;

    if (_flagEloss) {
      _hEloss = tfs->make<TTree>("hEloss", "Energy loss info");
      _hEloss->Branch("mean", &_eloss_mean, "mean/F");
      _hEloss->Branch("mp", &_eloss_mp, "mp/F");
      _hEloss->Branch("volid", &_eloss_vid, "volid/I");
      _hEloss->Branch("p", &_eloss_p, "p/F");
    }

    if (_flagDiagnostics) {
      _hExitCode = tfs->make<TH1F>("hExitCode", "Exit code", 10, -1, 9);
      _hNSteps = tfs->make<TH1F>("hNSteps", "Number of extrapolation steps", 300, 0, 30000./_extrapolationStep);
      _hNData = tfs->make<TH1F>("hNData", "Number of data points", 200, 0, _maxNBack);
      _hFL = tfs->make<TH1F>("hFL", "Flight length", 300, 0, 30.);
      _hNPAClust = tfs->make<TH1F>("hNPAClust", "Number of PA cluster", 20, 0, 20);
      _hNSTClust = tfs->make<TH1F>("hNSTClust", "Number of ST cluster", 20, 0, 20);
      _hPinit = tfs->make<TH1F>("hPinit", "Initial momentum", 200, 90, 110);
      _hPfin = tfs->make<TH1F>("hPfin", "Final momentum", 200, 90, 110);
    }

  }

  void TrkExt::beginRun(art::Run & run){
    cerr << "TrkExt: From beginRun: " << run.id().run() << endl;
    GeomHandle<DetectorSystem> det;
    _origin = det->toMu2e( CLHEP::Hep3Vector(0.,0.,0.) );   // add this to transfer detector coord. to mu2e coord.
    GeomHandle<WorldG4> g4world;
    _mu2eOriginInWorld = g4world->mu2eOriginInWorld(); // add this to transfer mu2e coord. to g4 coord. 
    cerr << "TrkExt: Detector coord origin in mu2e coord = " << _origin << endl;
    cerr << "TrkExt: Mu2e coord origin in G4 coord = " << _mu2eOriginInWorld << endl;


  }

  void TrkExt::beginSubRun(art::SubRun & lblock ) {
    cerr << "TrkExt: From beginSubRun. " << endl;
    _bfMgr = GeomHandle<BFieldManager>().get();
    _mydet.initialize();
  }

  void TrkExt::endJob(){
    cerr << "TrkExt: From endJob. " << endl;
  }


////////// Produce ///////////

  void TrkExt::produce(art::Event& event) {
    _evtid = event.id().event();
    if (_evtid%100 == 0) cerr << "TrkExt: event " << _evtid << endl;

    auto_ptr<TrkExtTrajCollection> trajcol(new TrkExtTrajCollection);

    art::Handle<KalRepCollection> trksHandle;
    event.getByLabel(_fitterModuleLabel, _trkPatRecInstanceName, trksHandle);
    KalRepCollection const& trks = *trksHandle;


    for ( size_t i=0; i< trks.size(); ++i ){
      _trkid = i;
      KalRep const& trk   = *(trks.at(i));
      Hep3Vector xstart, pstart, xstop, pstop;
      HepMatrix covstart(6,6,0);
      HepMatrix covstop(6,6,0);
      readTrkPatRec (trk, &xstart, &pstart, &xstop, &pstop, &covstart, &covstop);

      //cerr << "Track extrapolation at " << _evtid << ", track " << _trkid << endl;

      _traj.clear();
      if (_useVirtualDetector) {
        TrkHotList const* hits  = trk.hotList();
        if (!(readVD(event, hits))) {
          cerr << "TrkExt Warning: Cannot read VD at evt " << _evtid << ", trk " << i << ". Skipping" << endl;
          trajcol->push_back(_traj);
          continue;
        }
        if (_vdx[2] < -99998 || _vdy[2] < -99998 || _vdz[2] <-99998) {
          cerr << "TrkExt Warning: VD2 info not found at evt " << _evtid << ", trk " << i << ". Skipping" << endl;
          trajcol->push_back(_traj);
          continue;
        }
        xstart.set (_vdx[2], _vdy[2], _vdz[2]);
        pstart.set (_vdpx[2], _vdpy[2], _vdpz[2]);
      }

      int nsteps = doExtrapolation (xstart, pstart, covstart, false);
      if (_flagDiagnostics) {
        _hExitCode->Fill(_traj.exitCode());
        _hNSteps->Fill(nsteps);
        _hNData->Fill(_traj.size());
        _hFL->Fill(_traj.flightLength()*0.001);
        _hNPAClust->Fill(_traj.paHitIndex().size());
        _hNSTClust->Fill(_traj.stHitIndex().size());
        _hPinit->Fill(_traj.front().momentum().mag());
        _hPfin->Fill(_traj.back().momentum().mag());
      }

      trajcol->push_back(_traj);

    }  // end of trks loop

    event.put(trajcol);
  }


////////// Utility functions ///////////

  bool TrkExt::checkOutofReflectionLimit(Hep3Vector & x, Hep3Vector & p) {
    if(p.z() >=0) return false;
    Hep3Vector xx = x+_origin;
    if(xx.z() > _dummyStoppingTarget_z0+_dummyStoppingTarget_halfLength) return true;
    return false;
  }



////////// BField functions ///////////

  Hep3Vector TrkExt::getBField (Hep3Vector& x) {
    // x in Detector coordinate
    //mu2e::GeomHandle<BFieldManager> bfMgr;
    Hep3Vector b = _bfMgr->getBField(x + _origin);
    if (b.mag() >10) {
      Hep3Vector xx = x+_origin;
      cerr << "TrkExt: Crazy bfield : (" << b.x() << ", " << b.y() << ", " << b.z() << ") at (" << xx.x() << ", " << xx.y() << ", " << xx.z() << ")" << endl;
    }
    return b;
  }

  Hep3Vector TrkExt::getBField (HepVector& r) {
    Hep3Vector x(r[0], r[1], r[2]);
    return getBField(x);
  }


  Hep3Vector TrkExt::getBFieldWithGradient( Hep3Vector & x, 
                                      double & bxx, double & bxy, double & bxz, 
                                      double & byx, double & byy, double & byz, 
                                      double & bzx, double & bzy, double & bzz) {

    Hep3Vector B0 = getBField(x);

    if (_bFieldGradientMode == 1) {
      double xx = x.x();
      double yy = x.y();
      double zz = x.z();
      double h = 5.;
      Hep3Vector mx (xx-h,yy,zz);
      Hep3Vector px (xx+h,yy,zz);
      Hep3Vector my (xx,yy-h,zz);
      Hep3Vector py (xx,yy+h,zz);
      Hep3Vector mz (xx,yy,zz-h);
      Hep3Vector pz (xx,yy,zz+h);

      Hep3Vector Bmx = getBField(mx);
      Hep3Vector Bpx = getBField(px);
      Hep3Vector Bmy = getBField(my);
      Hep3Vector Bpy = getBField(py);
      Hep3Vector Bmz = getBField(mz);
      Hep3Vector Bpz = getBField(pz);

      bxx = (Bpx.x() - Bmx.x()) / (2.*h);
      bxy = (Bpy.x() - Bmy.x()) / (2.*h);
      bxz = (Bpz.x() - Bmz.x()) / (2.*h);
      byx = (Bpx.y() - Bmx.y()) / (2.*h);
      byy = (Bpy.y() - Bmy.y()) / (2.*h);
      byz = (Bpz.y() - Bmz.y()) / (2.*h);
      bzx = (Bpx.z() - Bmx.z()) / (2.*h);
      bzy = (Bpy.z() - Bmy.z()) / (2.*h);
      bzz = (Bpz.z() - Bmz.z()) / (2.*h);
    }
    else {
      bxx = 0;
      bxy = 0;
      bxz = 0;
      byx = 0;
      byy = 0;
      byz = 0;
      bzx = 0;
      bzy = 0;
      bzz = 0;
    }
    return B0;
  }


/////////// Read VD //////////////

  bool TrkExt::readVD (const art::Event& event, TrkHotList const * hits) {
    unsigned int i = 0, j, k = -1;

    // iterate from hot list
    for (TrkHotList::hot_iterator iter = hits->begin() ; iter != hits->end() ; ++iter) {
      ++k;
      const TrkHitOnTrk * hit  = iter.get();
      // read assoc. TrkStrawHit
      const mu2e::TrkStrawHit* trkStrawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hit);

      // try next if not found
      if (!trkStrawHit) continue;

      // read assoc. StrawHit
      const StrawHit& strawHit = trkStrawHit->strawHit();

      // read StrawHitCollection
      art::Handle<StrawHitCollection> shcHandle;
      event.getByLabel(_makerModuleLabel, shcHandle);
      StrawHitCollection const & shc = *shcHandle;

      // find the same StrawHit entry from StrawHitCollection
      int hitid = -1;
      for (i = 0 ; i <shc.size() ; ++i) {
        StrawHit  const& sh(shc.at(i));
        if (sh == strawHit) {
          hitid = int(i);
          break;
        }
      }

      // if not found, try next
      if (hitid <0) continue;

      // find StepPointMC assoc. to StrawHit
      art::Handle<PtrStepPointMCVectorCollection> stepsCHandle;
      event.getByLabel(_makerModuleLabel, "StrawHitMCPtr", stepsCHandle);
      PtrStepPointMCVectorCollection const & stepsC = *stepsCHandle;

      if (stepsC.size() <=0) continue;

      PtrStepPointMCVector const & steps (stepsC.at(hitid));

      if (steps.size() <=0) continue;

      StepPointMC const & step = *steps.at(0);

      // read SimParticle from StepPointMC
      SimParticle const &sim = *step.simParticle();

      // now simid is found    
      cet::map_vector_key simid =  sim.id();

      // read VD
      vector<Hep3Vector> vdp[5];
      vector<Hep3Vector> vdx[5];
      for (i = 0 ; i <5 ; ++i) {
        vdp[i].clear();
        vdx[i].clear();
      }
      art::Handle<StepPointMCCollection> hits;
      event.getByLabel(_g4ModuleLabel,"virtualdetector",hits);
      for ( i=0; i<hits->size(); ++i ){
        const StepPointMC& hit = (*hits)[i];
        if (hit.trackId() != simid) continue;
  
        switch (hit.volumeId()) {
          case VirtualDetectorId::ST_In :
            j = 0;
            break;
          case VirtualDetectorId::ST_Out :
            j = 1;
            break;
          case VirtualDetectorId::TT_FrontHollow :
          case VirtualDetectorId::TT_FrontPA :
            j = 2;
            break;
          case VirtualDetectorId::TT_MidInner :
          case VirtualDetectorId::TT_Mid :
            j = 3;
            break;
          case VirtualDetectorId::TT_Back :
            j = 4;
            break;
          default :
            j = 9999;
            break;
        }
        if (j >10) continue;
        vdp[j].push_back(hit.momentum());
        vdx[j].push_back(hit.position()-_origin);
      }
      for (i = 0 ; i <5 ; ++i) {
        if (vdp[i].size() >0) {
          Hep3Vector psum(0,0,0);
          Hep3Vector xsum(0,0,0);
          for (unsigned int k = 0 ; k < vdp[i].size() ; ++k) {
            psum += (vdp[i])[k];
            xsum += (vdx[i])[k];
          }
  
          psum /=  double(vdp[i].size());
          xsum /=  double(vdx[i].size());
          _vdx[i] = xsum.x();
          _vdy[i] = xsum.y();
          _vdz[i] = xsum.z();
          _vdpx[i] = psum.x();
          _vdpy[i] = psum.y();
          _vdpz[i] = psum.z();
          _vdp[i] = psum.mag();
        }
        else {
          _vdx[i] = -99999;
          _vdy[i] = -99999;
          _vdz[i] = -99999;
          _vdpx[i] = -99999;
          _vdpy[i] = -99999;
          _vdpz[i] = -99999;
          _vdp[i] = -99999;
        }
      }
      // read successful. return true
      return true;
    }

    // Reaching here means VD reading not successful. return false

    return false;
  }


/////////// Read TrkPatRec //////////////

  void TrkExt::readTrkPatRec(KalRep const & krep, Hep3Vector * xstart, Hep3Vector * pstart, Hep3Vector * xstop, Hep3Vector * pstop, HepMatrix * covstart, HepMatrix * covstop) {

    double _trkl0 = krep.startValidRange();
    double _trkl1 = krep.endValidRange();
    HepPoint xstart_ = krep.position(_trkl0);
    Hep3Vector pstart_ = krep.momentum(_trkl0);
    HepPoint xstop_ = krep.position(_trkl1);
    Hep3Vector pstop_ = krep.momentum(_trkl1);
    HepSymMatrix xxcov0(3,0);
    HepSymMatrix xxcov1(3,0);
    HepSymMatrix ppcov0(3,0);
    HepSymMatrix ppcov1(3,0);
    HepMatrix xpcov0(3,3,0);
    HepMatrix xpcov1(3,3,0);
    HepMatrix covstart_(6,6,0);
    HepMatrix covstop_(6,6,0);
    krep.getAllCovs(_trkl0, xxcov0, ppcov0, xpcov0);
    krep.getAllCovs(_trkl1, xxcov1, ppcov1, xpcov1);

    covstart_[0][0] = xxcov0[0][0];
    covstart_[0][0] = xxcov0[0][0];
    covstart_[0][1] = xxcov0[0][1];
    covstart_[0][2] = xxcov0[0][2];
    covstart_[0][3] = xpcov0[0][0];
    covstart_[0][4] = xpcov0[0][1];
    covstart_[0][5] = xpcov0[0][2];
    covstart_[1][1] = xxcov0[1][1];
    covstart_[1][2] = xxcov0[1][2];
    covstart_[1][3] = xpcov0[1][0];
    covstart_[1][4] = xpcov0[1][1];
    covstart_[1][5] = xpcov0[1][2];
    covstart_[2][2] = xxcov0[2][2];
    covstart_[2][3] = xpcov0[2][0];
    covstart_[2][4] = xpcov0[2][1];
    covstart_[2][5] = xpcov0[2][2];
    covstart_[3][3] = ppcov0[0][0];
    covstart_[3][4] = ppcov0[0][1];
    covstart_[3][5] = ppcov0[0][2];
    covstart_[4][4] = ppcov0[1][1];
    covstart_[4][5] = ppcov0[1][2];
    covstart_[5][5] = ppcov0[2][2];

    covstop_[0][0] = xxcov1[0][0];
    covstop_[0][1] = xxcov1[0][1];
    covstop_[0][2] = xxcov1[0][2];
    covstop_[0][3] = xpcov1[0][0];
    covstop_[0][4] = xpcov1[0][1];
    covstop_[0][5] = xpcov1[0][2];
    covstop_[1][1] = xxcov1[1][1];
    covstop_[1][2] = xxcov1[1][2];
    covstop_[1][3] = xpcov1[1][0];
    covstop_[1][4] = xpcov1[1][1];
    covstop_[1][5] = xpcov1[1][2];
    covstop_[2][2] = xxcov1[2][2];
    covstop_[2][3] = xpcov1[2][0];
    covstop_[2][4] = xpcov1[2][1];
    covstop_[2][5] = xpcov1[2][2];
    covstop_[3][3] = ppcov1[0][0];
    covstop_[3][4] = ppcov1[0][1];
    covstop_[3][5] = ppcov1[0][2];
    covstop_[4][4] = ppcov1[1][1];
    covstop_[4][5] = ppcov1[1][2];
    covstop_[5][5] = ppcov1[2][2];
    for (int i = 1 ; i <6 ; ++i) {
      for (int j = 0 ; j <i ; ++j) {
        covstart_[i][j] = covstart_[j][i];
        covstop_[i][j] = covstop_[j][i];
      }
    }

    xstart->set(xstart_.x(), xstart_.y(), xstart_.z()); 
    pstart->set(pstart_.x(), pstart_.y(), pstart_.z()); 
    xstop->set(xstop_.x(), xstop_.y(), xstop_.z()); 
    pstop->set(pstop_.x(), pstop_.y(), pstop_.z()); 
    *covstart = covstart_;
    *covstop  = covstop_;
    //see TrkMomCalculator::calcCurvAllCovs at BaBar/TrkBase/src/TrkMomCalculator.cc 
    return;
  }

//
///////// Track extrapolation ////////////

  int TrkExt::doExtrapolation (Hep3Vector xx, Hep3Vector pp, HepMatrix ccov, bool direction) {
    if (ccov.num_row() != 6 || ccov.num_col() != 6) {
      cerr << "TrkExt Warning : cannot use cov matrix." <<endl;
    }

    double stepSign = 1.;
    if (direction) stepSign = 1.;
    else stepSign = -1.;

    _extrapolationStep = stepSign * fabs(_extrapolationStep);

    double ds = _extrapolationStep;

    double s = 0;
    double dds = 0;
    double de, de_mean, de_mp, deltapp;

    TrkExtDetectorList::Enum prevolumeid = _mydet.volumeId(xx);
    TrkExtExitCode::Enum exitcode(TrkExtExitCode::Undefined);


    // Initialize TrajPoints
    //HepMatrix cov_init(6,6, 0);
    TrkExtTrajPoint r0 (0, xx, pp, ccov, TrkExtDetectorList::Enum(prevolumeid), 0, 0); // initial data
    TrkExtTrajPoint r1; // end data

    // Extrapolation stepping start
    int nsteps;
    for (nsteps = 0 ; ; ++nsteps) {

      // initial step size
      ds = _extrapolationStep;

      // Estimate next position
      r1 = calculateNextPosition(r0, ds);

      // check the volume info
      if (r1.volumeId() != TrkExtDetectorList::Undefined) {

        // if not same volume
        if (r1.volumeId() != r0.volumeId()) {

          // find intersection along the line of r0 and r1. The volume of intersection is always equal to that of r1
          Hep3Vector intersection = _mydet.intersection(r0.position(), r1.position());
      
          // new step size: to the intersection 
          ds = stepSign * ((intersection-r0.position()).mag());

          // estimate next position again
          r1 = calculateNextPosition(r0, ds);
          /*if (r1.volumeId() == r0.volumeId()) {
            do {
              ds += (stepSign * _mydet.limit());
              cout << "  small increment in ds " << endl;
              r1 = calculateNextPosition(r0, ds);
            } while (r1.volumeId() == r0.volumeId());
          }*/
          //cout << "  final ds = " << ds << endl;
        }
      } //end of volume check

      // flight length: global (s)  and local (dds)
      s += ds;
      dds += ds;

      // momentum loss term 
      deltapp = 0;

      // material effect check
      if (_turnOnMaterialEffect) {

        // get mean energy loss
        de_mean = _mydet.meanEnergyLoss(r0.momentum(), fabs(ds), TrkExtDetectorList::Enum(r0.volumeId()));

        // get most probable energy loss
        de_mp = _mydet.mostProbableEnergyLoss(r0.momentum(), fabs(ds), TrkExtDetectorList::Enum(r0.volumeId()));

        // deside which will be used
        if (_useStoppingPower) { de = de_mean; }
        else { de = de_mp; }

        // if energy loss do exist
        if (de >0) {

          // convert to momentum loss
          double e1 = safeSqrt(r0.momentum().mag()*r0.momentum().mag() + MEC22) - stepSign*de;
          double p1 = safeSqrt(e1*e1-MEC22);
          double sf = p1 / r0.momentum().mag();
          deltapp = sf -1.;

          // apply to the momentum of next position (r1)
          r1.scaleMomentum(sf);

          // record energy loss info
          if (_flagEloss) {
            _eloss_mean = de_mean / fabs(ds);
            _eloss_mp = de_mp / fabs(ds);
            _eloss_vid = r0.volumeId();
            _eloss_p = r0.momentum().mag();
            _hEloss->Fill();
          }
        }
      }  // end of material effect check

      // covariance calculation -calculating covariance from transport is default
      HepMatrix cov1 = getCovarianceTransport(r0, ds, deltapp);

      // calculate covariance from multiple scattering and add to previous one - it's optional
      if (_turnOnMultipleScattering) {
        HepMatrix cov2 = getCovarianceMultipleScattering(r0, ds);
        HepMatrix cov = cov1 + cov2;
        r1.setCovariance (cov);
      }
      else {
        r1.setCovariance (cov1);
      }

      // PA and ST hit booking
      if (r0.volumeId() != r1.volumeId()) {
        if      (r1.volumeId() == TrkExtDetectorList::ProtonAbsorber) _traj.addPAHit(r1.trajPointId());
        else if (r1.volumeId() == TrkExtDetectorList::StoppingTarget) _traj.addSTHit(r1.trajPointId());
      }

      // loop exit conditions : record and exit loop 
      if      ( r1.volumeId() == TrkExtDetectorList::Undefined )          exitcode = TrkExtExitCode::UndefinedVolume;
      else if ( r1.momentum().mag() > _maxMomentum )                      exitcode = TrkExtExitCode::MaximumMomentum;
      else if ( checkOutofReflectionLimit(r1.position(), r1.momentum()) ) exitcode = TrkExtExitCode::ReflectionLimit;
      else if ( int(_traj.size()) >= _maxNBack-1)                         exitcode = TrkExtExitCode::MaximumPoints;
      else if (   fabs(dds) > _recordingStep 
               || r0.volumeId() != r1.volumeId()
               || r0.volumeId() == TrkExtDetectorList::ProtonAbsorber
               || r0.volumeId() == TrkExtDetectorList::StoppingTarget   ) exitcode = TrkExtExitCode::WriteData;
      else                                                                exitcode = TrkExtExitCode::DoNotExit;

      if      (exitcode == TrkExtExitCode::DoNotExit) { ; }
      else if (exitcode == TrkExtExitCode::Undefined) { ; }
      else if (exitcode == TrkExtExitCode::WriteData) {
          _traj.push_back(r1);
          dds = 0;
      }
      else {
        _traj.push_back(r1);
        break;
      }
       
      // update start point
      r0 = r1;

    } // end of stepping

    // book track-wide variable 
    _traj.setExitCode(exitcode);
    return nsteps;
  }


///////// Covariance ////////////

  HepMatrix TrkExt::getCovarianceTransport(TrkExtTrajPoint & r0, double ds, double deltapp) {
    HepMatrix & E = r0.covariance();
    HepMatrix Ep(6,6,0);
    if (E.num_row() !=6 || E.num_col() !=6) {
      cerr << "TrkExt Warning : cannot calculate covariance" << endl;
      return Ep;
    }
    HepMatrix J(6,6,0);
    double px = r0.px();
    double py = r0.py();
    double pz = r0.pz();
    double p = r0.momentum().mag();
    if (p == 0) {
      cerr << "TrkExt Warning : 0 momentum?" << endl;
      return Ep;
    }
    double pp = p*p;
    double ppp = pp*p;

    double Bx;
    double By;
    double Bz;
    double Bxx;
    double Bxy;
    double Bxz;
    double Byx;
    double Byy;
    double Byz;
    double Bzx;
    double Bzy;
    double Bzz;

    Hep3Vector B = getBFieldWithGradient (r0.position(), Bxx, Bxy, Bxz, Byx, Byy, Byz, Bzx, Bzy, Bzz) ;
    Bx = B.x();
    By = B.y();
    Bz = B.z();

    double kqds = ds * RUNGE_KUTTA_KQ;

    J[0][0] = 1;
    J[0][1] = 0;
    J[0][2] = 0;
    J[0][3] = ds*(py*py+pz*pz)/ppp;
    J[0][4] = -ds*px*py/ppp;
    J[0][5] = -ds*px*pz/ppp;

    J[1][0] = 0;
    J[1][1] = 1;
    J[1][2] = 0;
    J[1][3] = -ds*py*px/ppp;
    J[1][4] = ds*(pz*pz+px*px)/ppp;
    J[1][5] = -ds*py*pz/ppp;

    J[2][0] = 0;
    J[2][1] = 0;
    J[2][2] = 1;
    J[2][3] = -ds*pz*px;
    J[2][4] = -ds*pz*py;
    J[2][5] = ds*(px*px+py*py)/ppp;

    J[3][0] = kqds /p *(py*Bzx - pz*Byx);
    J[3][1] = kqds /p *(py*Bzy - pz*Byy);
    J[3][2] = kqds /p *(py*Bzz - pz*Byz);
    J[3][3] = 1+deltapp-kqds/ppp*px*(py*Bz-pz*By);
    J[3][4] = kqds/ppp *( Bz*pp - py*(py*Bz-pz*By));
    J[3][5] = kqds/ppp *(-By*pp - pz*(py*Bz-pz*By));

    J[4][0] = kqds /p *(pz*Bxx - px*Bzx);
    J[4][1] = kqds /p *(pz*Bxy - px*Bzy);
    J[4][2] = kqds /p *(pz*Bxz - px*Bzz);
    J[4][3] = kqds/ppp *(-Bz*pp - px*(pz*Bx-px*Bz));
    J[4][4] = 1+deltapp-kqds/ppp*py*(pz*Bx-px*Bz);
    J[4][5] = kqds/ppp *( Bx*pp - pz*(pz*Bx-px*Bz));

    J[5][0] = kqds /p *(px*Byx - py*Bxx);
    J[5][1] = kqds /p *(px*Byy - py*Bxy);
    J[5][2] = kqds /p *(px*Byz - py*Bxz);
    J[5][3] = kqds/ppp *( By*pp - px*(px*By-py*Bx));
    J[5][4] = kqds/ppp *(-Bx*pp - py*(px*By-py*Bx));
    J[5][5] = 1+deltapp-kqds/ppp*pz*(px*By-py*Bx);

    HepMatrix JT = J.T();

    Ep =  J*E*JT;

    return Ep;
  }

  HepMatrix TrkExt::getCovarianceMultipleScattering(TrkExtTrajPoint & r0, double ds) {
    Hep3Vector e = r0.momentum().unit();
    double costh1 = e.x();
    double costh2 = e.y();
    double costh3 = e.z();
    double sinth12  = 1. - costh1*costh1;
    double sinth22  = 1. - costh2*costh2;
    double sinth32  = 1. - costh3*costh3;

    HepMatrix R(3,3,0);
    R[0][0] = sinth12;
    R[0][1] = -costh1*costh2;
    R[0][2] = -costh1*costh3;
    R[1][0] = R[0][1];
    R[1][1] = sinth22;
    R[1][2] = -costh2*costh3;
    R[2][0] = R[0][2];
    R[2][1] = R[1][2];
    R[2][2] = sinth32;

    double th = _mydet.scatteringAngle(r0.momentum(), ds, TrkExtDetectorList::Enum(r0.volumeId()));
    double p = r0.momentum().mag();
    double m11 = 0.3333*th*th*ds*ds;
    double m12 = 0.5*th*th*ds*p;
    double m22 = th*th*p*p;
    
    HepMatrix Em(6,6,0);

    for (int i = 0 ; i < 3 ; ++i) {
      for (int j = 0 ; j <3 ; ++j) {
        Em[i][j] = m11 * R[i][j];
        Em[i+3][j] = Em[i][j+3] = m12 * R[i][j];
        Em[i+3][j+3] = m22 * R[i][j];
      }
    }

    return Em;
  }





///////// Functions for Runge-Kutta method ////////////

  TrkExtTrajPoint TrkExt::calculateNextPosition (TrkExtTrajPoint r00, double ds) {
    HepVector r0(6), re(6);
    r0 = r00.vector();
    HepVector dr1_ds = _runge_kutta_newpar_f(r0, getBField(r0));   HepVector r1 = r0+0.5*ds*dr1_ds;
    HepVector dr2_ds = _runge_kutta_newpar_f(r1, getBField(r1));   HepVector r2 = r0+0.5*ds*dr2_ds;
    HepVector dr3_ds = _runge_kutta_newpar_f(r2, getBField(r2));   HepVector r3 = r0+ds*dr3_ds;
    HepVector dr4_ds = _runge_kutta_newpar_f(r3, getBField(r3));

    re = r0 + (dr1_ds/6. + dr2_ds/3. + dr3_ds/3. + dr4_ds/6.)*ds;

    Hep3Vector x(re[0], re[1], re[2]);
    int volid = _mydet.volumeId(x);

    double p = r00.momentum().mag();
    double v = p/safeSqrt(p*p+MEC22)*VELOCITY_OF_LIGHT;
    double ft = fabs(ds) / v * 1.e6;

    return TrkExtTrajPoint(r00.trajPointId()+1, re, volid, r00.flightLength()+fabs(ds), ft);
  }



  HepVector TrkExt::_runge_kutta_newpar_5th (HepVector r0, double ds, bool mode) {

//    static double a2 = 0.2;
//    static double a3 = 0.3;
//    static double a4 = 0.6;
//    static double a5 = 1.;
//    static double a6 = 0.825;
    static double b21 = 0.2;
    static double b31 = 0.075;
    static double b32 = 0.225;
    static double b41 = 0.3;
    static double b42 = -0.9;
    static double b43 = 1.2;
    static double b51 = -11./54.;
    static double b52 = 2.5;
    static double b53 = -70./27.;
    static double b54 = 35./27.;
    static double b61 = 1631./55296.;
    static double b62 = 175./512.;
    static double b63 = 575./13824.;
    static double b64 = 44275./110592.;
    static double b65 = 253./4096.;
    static double c1 = 37./378.;
    static double c2 = 0.;
    static double c3 = 250./621.;
    static double c4 = 125./594.;
    static double c5 = 0.;
    static double c6 = 512./1771.;
    static double c1s = 2825./27648.;
    static double c2s = 0.;
    static double c3s = 18575./48384.;
    static double c4s = 13525./55296.;
    static double c5s = 277./14336.;
    static double c6s = 0.25;

    HepVector k1 = ds*_runge_kutta_newpar_f(r0, getBField(r0)); HepVector r1 = r0 + b21*k1;
    HepVector k2 = ds*_runge_kutta_newpar_f(r1, getBField(r1)); HepVector r2 = r0 + b31*k1 + b32*k2;
    HepVector k3 = ds*_runge_kutta_newpar_f(r2, getBField(r2)); HepVector r3 = r0 + b41*k1 + b42*k2 + b43*k3;
    HepVector k4 = ds*_runge_kutta_newpar_f(r3, getBField(r3)); HepVector r4 = r0 + b51*k1 + b52*k2 + b53*k3 + b54*k4;
    HepVector k5 = ds*_runge_kutta_newpar_f(r4, getBField(r4)); HepVector r5 = r0 + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5;
    HepVector k6 = ds*_runge_kutta_newpar_f(r5, getBField(r5)); 

    if (mode) {
      return r0 + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6;
    }else {
      return r0 + c1s*k1 + c2s*k2 + c3s*k3 + c4s*k4 + c5s*k5 + c6s*k6;
    }
  }


  HepVector TrkExt::_runge_kutta_newpar_f (HepVector r, Hep3Vector B) {
    HepVector ret(6);
    Hep3Vector p(r[3], r[4], r[5]);
    Hep3Vector e = p.unit();
    Hep3Vector p_ = RUNGE_KUTTA_KQ* (e.cross(B));
    ret[0] = e.x();
    ret[1] = e.y();
    ret[2] = e.z();
    ret[3] = p_.x();
    ret[4] = p_.y();
    ret[5] = p_.z();
    return ret;
  }



} // end namespace mu2e

using mu2e::TrkExt;
DEFINE_ART_MODULE(TrkExt);
