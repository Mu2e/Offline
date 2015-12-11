//
// $Id: TrkExtDiag.hh,v 1.3 2013/05/16 18:23:39 mjlee Exp $
// $Author: mjlee $ 
// $Date: 2013/05/16 18:23:39 $
//
// Class for reading TrkExt
//
// Note on index for (EXT)VDHitInfo struct : 
//
//   0 for down-going (target->tracker) (pz>=0),
//   1 for up-going (tracker->target) (pz<0).
//
// See the note at TrkExtInstanceName.hh
// Note that direction of extrapolation stepping is opposite 
// for updown=0 track (ex. normal ce without reflection). 
// Don't be confused with extrapolation stepping direction.   
//
// Then, status is defined by :
//   0 : no data avaliable
//   1 : down going only
//   2 : up going only
//   3 : both avaliable
//
// Note that, for EXTVDHitInfo struct with normal ce extrapolation, 
// status = 1 is possible case (when there is no reflection), 
// when status = 2 is not possible.
// Mostly it should be status = 3.  
//
#ifndef TrkExtDiag_HH
#define TrkExtDiag_HH

#include "art/Framework/Principal/fwd.h"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "TTree.h"
#include "Rtypes.h"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{  
  const unsigned int MAXNHOT = 100;
  const unsigned int MAXNTRK = 100;
  const unsigned int MAXNSIM = 1000;
  const unsigned int MAXNEXT = 10000;
  const unsigned int MAXNPA = 50;
  const unsigned int MAXNST = 50;
  const unsigned int MAXNVD = 10;

  struct PositionInfo {
    Float_t x;
    Float_t y;
    Float_t z;
    PositionInfo() : x(0), y(0), z(0) {}
  };

  struct HitInfo {
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t t;
    Float_t s;
    HitInfo() : s(0) {}
  };

  struct VDHitInfo {
    Float_t x[2];
    Float_t y[2];
    Float_t z[2];
    Float_t px[2];
    Float_t py[2];
    Float_t pz[2];
    Float_t p[2];
    Float_t d0[2];
    Float_t p0[2];
    Float_t om[2];
    Float_t z0[2];
    Float_t td[2];
    UInt_t status;
    VDHitInfo() : status(0) {}
    void clear() {
      for (int i = 0 ; i <2 ; ++i) {
        x[i] = y[i] = z[i] = px[i] = py[i] = pz[i] = 0;
        p[i] = d0[i] = p0[i] = om[i] = z0[i] = td[i] = 0;
        status = 0;
      }
    }
  };

  struct EXTVDHitInfo {
    Float_t r[2];
    Float_t p[2];
    Float_t d0[2];
    Float_t p0[2];
    Float_t om[2];
    Float_t z0[2];
    Float_t td[2];
    Float_t er[2];
    Float_t ep[2];
    Float_t ed0[2];
    Float_t ep0[2];
    Float_t eom[2];
    Float_t ez0[2];
    Float_t etd[2];
    UInt_t status;
    EXTVDHitInfo() : status(0) {}
    void clear() {
      for (int i = 0 ; i <2 ; ++i) {
        r[i] = p[i] = d0[i] = p0[i] = om[i] = z0[i] = td[i] = 0;
        er[i] = ep[i] = ed0[i] = ep0[i] = eom[i] = ez0[i] = etd[i] = 0;
        status = 0;
      }
    }
  };
  
  struct HelixParameterInfo {
    Float_t d0;
    Float_t p0;
    Float_t om;
    Float_t z0;
    Float_t td;
    Float_t ed0;
    Float_t ep0;
    Float_t eom;
    Float_t ez0;
    Float_t etd;
    HelixParameterInfo() : 
      d0(0.), p0(0.), om(0.), z0(0.), td(0.),
      ed0(0.), ep0(0.), eom(0.), ez0(0.), etd(0.) {}
    void clear() {
      d0 = p0 = om = z0 = td = 0;
      ed0 = ep0 = eom = ez0 = etd = 0;
    }
  };


  class TrkExtDiag {
  public:
    explicit TrkExtDiag(fhicl::ParameterSet const&);
    virtual ~TrkExtDiag();
    TTree* createTrkExtDiag();
    void setRunInfo();
    void setSubRunInfo();
    void trkExtDiag(const art::Event & evt, const KalRep & trk, const TrkExtTraj & trkext);
    void trkExtDiag(void);

  private:
    unsigned int readHit(const art::Event &evt, const TrkHitVector & hot);
    unsigned int readTrk(const art::Event &evt, const KalRep & krep);
    unsigned int readMC(const art::Event &evt, const KalRep & krep, const TrkHitVector & hot);
    unsigned int readExt(const art::Event &evt, const TrkExtTraj & trkext);

    CLHEP::HepVector getHelixParameters (const CLHEP::Hep3Vector & x, const CLHEP::Hep3Vector & p, int sign) const ;
    CLHEP::HepVector getHelixParametersErr (const CLHEP::Hep3Vector & x, const CLHEP::Hep3Vector & p,  CLHEP::Hep3Vector &ex, CLHEP::Hep3Vector &ep, int sign) const ;
    double findTurnAround (double s1, double s2, double s3, double f1, double f2, double f3);
    double findTurnAroundSim (int i, double f1, double f2, double f3);
    double interpolate (double s, double s1, double s2, double s3, double f1, double f2, double f3);
    double interpolate2 (double s, double s1, double s2, double f1, double f2);
    double getRadialError (const CLHEP::Hep3Vector & x, const CLHEP::HepMatrix & cov);
    double getMomentumError (const CLHEP::Hep3Vector & p, const CLHEP::HepMatrix & cov);

  private:
    // event-wide 
    int _evtid, _hepid, _updown, _trkid, _exitcode;
    // hit
    unsigned int _nhots;
    float _hotx[MAXNHOT], _hoty[MAXNHOT], _hotz[MAXNHOT], _hott0[MAXNHOT];
    // trk
    unsigned int _ntrks;
    float _trkl0, _trkl1, _trkl[MAXNTRK], _trkx[MAXNTRK], _trky[MAXNTRK], _trkz[MAXNTRK];
    float _trkpx[MAXNTRK], _trkpy[MAXNTRK], _trkpz[MAXNTRK], _trkp[MAXNTRK];
    HelixParameterInfo _trk0;
    HelixParameterInfo _trk1;
    // mc
    unsigned int _nsim;
    float _simx[MAXNSIM], _simy[MAXNSIM], _simz[MAXNSIM], _simp0, _simt0;
    // mc - turn around point
    PositionInfo _simtp;
    int _simtpqual;
    // mc - pa and st
    unsigned int _nmcpa, _nmcst;
    float _mcpapx[MAXNPA], _mcpapy[MAXNPA], _mcpapz[MAXNPA], _mcpap[MAXNPA], _mcpadp[MAXNPA], _mcpade[MAXNPA], _mcpadei[MAXNPA], _mcpadeni[MAXNPA], _mcpaz[MAXNPA];
    float _mcstpx[MAXNST], _mcstpy[MAXNST], _mcstpz[MAXNST], _mcstp[MAXNST], _mcstdp[MAXNST], _mcstde[MAXNST], _mcstdei[MAXNST], _mcstdeni[MAXNST], _mcstx[MAXNPA], _mcsty[MAXNPA], _mcstz[MAXNPA], _mcstt[MAXNPA];
    // mc - vd
    VDHitInfo _vdsi;
    VDHitInfo _vdso;
    VDHitInfo _vdtf;
    VDHitInfo _vdtm;
    VDHitInfo _vdtb;
    // trkext 
    unsigned int _next;
    float _extx[MAXNEXT], _exty[MAXNEXT], _extz[MAXNEXT];
    float _extpx[MAXNEXT], _extpy[MAXNEXT], _extpz[MAXNEXT], _extp[MAXNEXT];
    float _extrho[MAXNEXT], _exts[MAXNEXT], _extt[MAXNEXT];
    float _extex[MAXNEXT], _extey[MAXNEXT], _extez[MAXNEXT];
    float _extepx[MAXNEXT], _extepy[MAXNEXT], _extepz[MAXNEXT], _extep[MAXNEXT], _exter[MAXNEXT]; 
    int _extvid[MAXNEXT];
    // trkext - turn around point
    int _extitp;
    HitInfo _exttp;
    // trkext - pa and st
    unsigned int _nextpa, _nextst;
    float _extpaz[MAXNPA], _extpadp[MAXNPA], _extpadptot;
    float _extstx[MAXNST], _extsty[MAXNST], _extstz[MAXNST], _extstt[MAXNST], _extstdp[MAXNST], _extstdptot;
    // trkext - vd
    EXTVDHitInfo _extvdsi;
    EXTVDHitInfo _extvdso;
    EXTVDHitInfo _extvdtf;
    EXTVDHitInfo _extvdtm;
    EXTVDHitInfo _extvdtb;

    //local variables 
    std::string _makerModuleLabel;
    std::string _g4ModuleLabel;
    CLHEP::Hep3Vector _origin;
    double _vdzsi, _vdzso, _vdztf, _vdztm, _vdztb;
    BFieldManager const * _bfMgr;


  public:
    TTree *_extdiag;
  };
}

#endif

