//
// output utilities for reco modules
//
// $Id: kalFitDataOuts.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//

//C++
#include <cstddef>

//ROOT
#include "Rtypes.h"

// data
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

namespace mu2e 
{

struct HitInfo {
  CLHEP::Hep3Vector _pos;
  Float_t _edep, _disttomid, _disttomid_err;//, _time, _corrtime;
  Int_t _Slayer, _layer, _cell;
  //Bool_t _vloose, _loose, _tight, _delta;
  CLHEP::Hep3Vector _mcpos;
  Int_t _mcpdg, _mcgen, _mcproc;
  Float_t _mcedep, _mcdisttomid;//, _mctime, _mct0, _mcmom, _mctd;
// root
  ClassDef(HitInfo,1)
};

// struct to describe an arc within a vector of TrkStrawHits
 struct TrkArc {
   size_t _begin;
   size_t _end;
   size_t _beginactive;
   size_t _endactive;
   unsigned _ntsh;
   unsigned _nactive;
// initialize to 0
   TrkArc(size_t begin=0) : _begin(begin),_end(begin),_beginactive(begin),_endactive(begin),_ntsh(0),_nactive(0){}
 };

 struct TrkArcInfo {
   Int_t _narctsh, _narcactive;
   Float_t _arctshlen, _arcactivelen;
// root
   ClassDef(TrkArcInfo,1)
 };

struct TrkCellHitInfo {
  Int_t _active, _usable, _Slayer, _layer, _cell;
  Float_t _z, _phi, _rho;
  Float_t _resid, _residerr, _rdrift, _rdrifterr, _trklen;
  Float_t _doca, _exerr, _penerr, _t0err;
  Float_t _ht, _tddist, _tdderr, _hlen;
  Int_t _ambig;
  Int_t _iarc, _iarchit;
  Float_t _architlen, _gaplow, _gaphi;
  Int_t _mcn, _mcnunique, _mcppdg, _mcpgen, _mcpproc;
  Int_t _mcpdg, _mcgen, _mcproc;
  Float_t _mcht, _mcdist, _mclen;
  Int_t _mcambig;
// root
  ClassDef(TrkCellHitInfo,1)
};

struct MCEvtData {
  MCEvtData(const StrawHitMCTruthCollection* mcstrawhits,
    const PtrStepPointMCVectorCollection* mchitptr,
    const StepPointMCCollection *mcsteps,
    const StepPointMCCollection *mcvdsteps) : _mcstrawhits(mcstrawhits),_mchitptr(mchitptr),
    _mcsteps(mcsteps),_mcvdsteps(mcvdsteps){}
  void clear() {_mcstrawhits = 0; _mchitptr = 0; _mcsteps = 0; _mcvdsteps = 0; _simparts = 0; }
  MCEvtData() {clear();}
  bool good() { return _mcstrawhits != 0 && _mchitptr != 0 && _mcsteps != 0 && _mcvdsteps != 0; }
  const StrawHitMCTruthCollection* _mcstrawhits;
  const PtrStepPointMCVectorCollection* _mchitptr;
  const StepPointMCCollection *_mcsteps, *_mcvdsteps;
  const SimParticleCollection *_simparts;
};

struct RecoInfo {
        int fit;
        float chi2,chi2in,chi2out,fitcon,radlen;
        int ndof,niter,ninter,nactive,nsites;
        int nhits,nhitstot,nturn;
        float genmom,wgt,momtrackerin,momin,momout;
        float t0;
        float fitmom, fitmomerr, seedmom;
        float fitmomout,fitmombeam;
        int fitnturn;
        float t0init,t0fit,errt0;
        int iseed;
        int iev;
        float trkfrstHitPosX, trkfrstHitPosY, trkfrstHitPosZ;
        float trk1Hitflght,trk1SimHitflght,trkto1SimHitdist;
        int atSimHit;
        int genID;
        int fromTargets;
        float lowFitRange, hiFitRange, firsthitPLength;
        float frstHitPosX, frstHitPosY, frstHitPosZ;
        int nHitSwire, nhitFwire;
        float pathInSwire, pathInFwire, pathInGas, simPathInGas, simPathInGasFrstLoop, simTotPath;
        float simMomLossInVol, simMomLossInVolFrstLoop, simMomLossInIWFrstloop;
        void clear(){
                fit=0; chi2=-1; chi2in=-1; chi2out=-1; fitcon=-1; radlen=-1;
                ndof=0; niter=0; ninter=0; nactive=0; nsites=0; nhits=0; nhitstot=0; nturn=0;
                genmom=0; wgt=1.0;momtrackerin=0; momin=0; momout=0; t0=0;
                fitmom=0; fitmomerr=0; seedmom=0; fitmomout=0;fitmombeam=0; fitnturn=0;
                t0init=0; t0fit=0; errt0=-1; iseed=-1; iev=-1;
                genID=-1; fromTargets=0; lowFitRange=0; hiFitRange=0; firsthitPLength=0; frstHitPosX=0; frstHitPosY=0; frstHitPosZ=0;
                trkfrstHitPosX=0; trkfrstHitPosY=0; trkfrstHitPosZ=0;
                trk1Hitflght=0;trk1SimHitflght=0;trkto1SimHitdist=-1;
                atSimHit=0;
                nHitSwire=0; nhitFwire=0; pathInSwire=0.0; pathInFwire=0.0;
                pathInGas=0.0; simPathInGas=0.0; simPathInGasFrstLoop=0.0; simTotPath=0.0;
                simMomLossInVol=0.0; simMomLossInVolFrstLoop=0.0; simMomLossInIWFrstloop=0.0;
        };
        void clearrec(){
                fit=0; chi2=-1; chi2in=-1; chi2out=-1; fitcon=-1; radlen=-1; ndof=0; niter=0; ninter=0; nactive=0; nsites=0;
                fitmom=0; fitmomerr=0; seedmom=0; fitmomout=0; fitmombeam=0;t0fit=0; t0init=0; errt0=-1; iseed=-1;
                trkfrstHitPosX=0; trkfrstHitPosY=0; trkfrstHitPosZ=0;
                trk1Hitflght=0;trk1SimHitflght=0;trkto1SimHitdist=-1;
                atSimHit=0;
        };

};

struct EventInfo {
        float elosstgt;
        int ntgt;
        void clear(){elosstgt=0;ntgt=0;};
};

}  // end namespace mu2e
