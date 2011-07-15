//
// MC functions associated with KalFit
// $Id: KalFitMC.hh,v 1.3 2011/07/15 04:44:06 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/07/15 04:44:06 $
//
#ifndef KalFitMC_HH
#define KalFitMC_HH

// data
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
#include "KalmanTests/inc/KalFit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "Rtypes.h"
#include "TTree.h"
#include "TClass.h"
#include <vector>

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

  struct TrkStrawHitInfo {
    Int_t _active,_usable;
    UInt_t _nmc;
// root macro
    ClassDef(TrkStrawHitInfo,1)
  };

  struct MCEvtData {
    MCEvtData(const StrawHitMCTruthCollection* mcstrawhits,
      const PtrStepPointMCVectorCollection* mchitptr,
      const StepPointMCCollection *mcsteps,
      const StepPointMCCollection *mcvdsteps) : _mcstrawhits(mcstrawhits),_mchitptr(mchitptr),
      _mcsteps(mcsteps),_mcvdsteps(mcvdsteps){}
    void clear() {_mcstrawhits = 0; _mchitptr = 0; _mcsteps = 0; _mcvdsteps = 0; }
    MCEvtData() {clear();}
    bool good() { return _mcstrawhits != 0 && _mchitptr != 0 && _mcsteps != 0 && _mcvdsteps != 0; }
    const StrawHitMCTruthCollection* _mcstrawhits;
    const PtrStepPointMCVectorCollection* _mchitptr;
    const StepPointMCCollection *_mcsteps, *_mcvdsteps;
  };
  
  typedef StepPointMCCollection::const_iterator MCStepItr;
//  struct test : public binary_function<double,double,bool> {
//    bool operator()(double x, double y) { return x < y;}
//  };

//  Simple helper class to find MC information within collections
  
  class KalFitMC {
  public:
    explicit KalFitMC(fhicl::ParameterSet const&);
    virtual ~KalFitMC();
// create a track definition object based on MC true particle
    bool trkFromMC(MCEvtData const& mcdata, cet::map_vector_key const& trkid, TrkDef& mytrk);
// diagnostic comparison of reconstructed tracks with MC truth
    void trkDiag(MCEvtData const& mcdata, TrkDef const& mytrk, TrkKalFit const& myfit);
    void hitDiag(MCEvtData const& mcdata, const TrkStrawHit* strawhit);
// allow creating the trees
    TTree* createTrkDiag();
    TTree* createHitDiag();
  private:
// helper functions
    void findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids,
      std::vector<MCStepItr>& steps);
// config parameters
    double _mintrkmom; // minimum true momentum at z=0 to create a track from
    double _mct0err;
    int _debug;
    unsigned _minnhits,_maxnhits;
    bool _purehits;
// vector of detector Ids corresponding to entrance and midplane
    std::vector<int> _midvids;
    std::vector<int> _entvids;
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
    std::vector<TrkStrawHitInfo> _tshinfo;

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
    Float_t _pdist;
    Float_t _pperp;
    Float_t _pmom;
  };
}

#endif

