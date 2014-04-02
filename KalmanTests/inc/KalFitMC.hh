//
// MC functions associated with KalFit
// $Id: KalFitMC.hh,v 1.36 2014/04/02 14:18:00 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/04/02 14:18:00 $
//
#ifndef KalFitMC_HH
#define KalFitMC_HH

// data
#include "art/Framework/Principal/fwd.h"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BaBar/BaBar.hh"
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
 // some convenient typedefs    
  typedef art::Ptr<SimParticle> SPPtr;
// structs
  struct MCHitSum {
    double _esum;
    unsigned _count;
    art::Ptr<SimParticle> _spp;
    int _pdgid;
    int _gid,_pid;
    StrawIndex _sid;
    double _t0, _time;
    CLHEP::Hep3Vector _pos;
    CLHEP::Hep3Vector _mom;
    MCHitSum() : _esum(0.0),_count(0),_pdgid(0),_gid(0), _sid(0), _time(0.0) {}
    MCHitSum(StepPointMC const& mchit,art::Ptr<SimParticle>& spp);
    MCHitSum(StrawDigiMC const& mcdigi);
    bool operator == (const MCHitSum& other) { return _spp == other._spp; }
    bool operator < (const MCHitSum& other) { return _spp < other._spp; }
    void append(StepPointMC const& mchit);
// comparison functor for ordering according to energy
    struct ecomp : public std::binary_function<MCHitSum,MCHitSum, bool> {
      bool operator()(MCHitSum const& t1, MCHitSum const& t2) { return t1._esum > t2._esum; }
    };
  };
  typedef std::vector<MCHitSum> MCHitSumVec;
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
    Int_t _active, _usable, _device, _sector, _layer, _straw;
    Float_t _z, _phi, _rho;
    Float_t _resid, _residerr, _rdrift, _rdrifterr, _trklen;
    Float_t _doca, _exerr, _penerr, _t0, _t0err;
    Float_t _ht, _tddist, _tdderr, _hlen;
    Float_t _edep, _dx;
    Int_t _ambig;
    Int_t _mcn, _mcnunique, _mcppdg, _mcpgen, _mcpproc;
    Int_t _mcpdg, _mcgen, _mcproc;
    Float_t _mct0, _mcht, _mcdist, _mclen;
    Float_t _mcedep;
    Int_t _mcambig;
    Bool_t _xtalk;
// root 
    ClassDef(TrkStrawHitInfo,1)
  };

  struct MCTrkInfo {
    Int_t _pdgid;
    Float_t _time;
    Float_t _mom;
    threevec _pos;
    helixpar _hpar;
    MCTrkInfo() : _pdgid(0), _time(0.0),_mom(0.0) {}
  };

  struct MCEvtData {
    MCEvtData(
      const PtrStepPointMCVectorCollection* mchitptr,
      const StepPointMCCollection *mcsteps,
      const StepPointMCCollection *mcvdsteps) : _mchitptr(mchitptr),
      _mcsteps(mcsteps),_mcvdsteps(mcvdsteps){}
    void clear() { _mchitptr = 0; _mcsteps = 0; _mcvdsteps = 0; _simparts = 0; _mcdigis = 0; }
    MCEvtData() {clear();}
    bool good() { return _mchitptr != 0 && _mcsteps != 0 && _mcvdsteps != 0; }
    const PtrStepPointMCVectorCollection* _mchitptr;
    const StepPointMCCollection *_mcsteps, *_mcvdsteps;
    const StrawDigiMCCollection *_mcdigis;
    const SimParticleCollection *_simparts;
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

  typedef StepPointMCCollection::const_iterator MCStepItr;
//  struct test : public binary_function<double,double,bool> {
//    bool operator()(double x, double y) { return x < y;}
//  };

//  Simple helper class to find MC information within collections
  
  class KalFitMC {
  public:
    enum TRACKERPOS {trackerEnt=0,trackerMid, trackerExit};
    enum relation {none=-1,same,daughter,mother,sibling,udaughter,umother,usibling};
    explicit KalFitMC(fhicl::ParameterSet const&);
    virtual ~KalFitMC();
// find MC data in the event.  This must be called each event, before the other functions
    bool findMCData(const art::Event& evt);
// create a track definition object based on MC true particle
    bool trkFromMC(cet::map_vector_key const& trkid, TrkDef& mytrk);
// diagnostic comparison of reconstructed tracks with MC truth
    void kalDiag(const KalRep* krep,bool fill=true);
    void hitsDiag(std::vector<const TrkStrawHit*> const& hits);
    void arcsDiag(std::vector<const TrkStrawHit*> const& hits);
    void mcTrkInfo(SimParticle const& sp);
    void hitDiag(const TrkStrawHit* strawhit);
// find associated sim particles to a track
    void findMCTrk(const KalRep* krep,art::Ptr<SimParticle>& spp);
// allow creating the trees
    TTree* createTrkDiag();
    TTree* createHitDiag();
// General StrawHit MC functions: these should move to more general class, FIXME!!
    const MCHitSumVec& mcHitSummary(size_t ihit) const { return _mchitsums[ihit]; }
// access to MC data
    MCEvtData const& mcData() const { return _mcdata; }
    void findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids,
	std::vector<MCStepItr>& steps);
// access to event-specific MC truth for conversion electron
    double MCT0(TRACKERPOS tpos) const;
    double MCMom(TRACKERPOS tpos) const;
    std::vector<int> const& VDids(TRACKERPOS tpos) const;
    const helixpar& MCHelix(TRACKERPOS tpos) const;
    double MCBrems() const { return _bremsesum; }
    static void fillMCHitSum(PtrStepPointMCVector const& mcptr,MCHitSumVec& summary );
// MC info about a track
    void fillMCTrkInfo(MCStepItr const& imcs, MCTrkInfo& trkinfo) const;
    void fillMCTrkInfo(SimParticle const& sp, MCTrkInfo& einfo) const;
    static relation relationship(art::Ptr<SimParticle> const& sppi,art::Ptr<SimParticle> const& sppj);
    static relation relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2);

  private:
// cache of event data
    MCEvtData _mcdata;
    const StrawHitCollection* _strawhits;
// cache of hit summary
   std::vector<MCHitSumVec> _mchitsums;
// event data labels
    std::string _mcstrawhitslabel;
    std::string _mcptrlabel;
    std::string _mcstepslabel;
    std::string _simpartslabel;
    std::string _mcdigislabel;
    std::string _strawhitslabel;
// helper functions
    static void findRelatives(PtrStepPointMCVector const& mcptr,std::map<SPPtr,SPPtr>& mdmap );
    void fillMCHitSummary();
    void fillHitsVector(const KalRep* krep,std::vector<const TrkStrawHit*>& hits);
    void findArcs(std::vector<const TrkStrawHit*> const& straws, std::vector<TrkArc>&  arcs) const;
    static int findArc(size_t itsh,std::vector<TrkArc>& arcs );
// config parameters
    double _mintrkmom; // minimum true momentum at z=0 to create a track from
    double _mct0err;
    bool _mcambig;
    int _debug,_diag;
    unsigned _minnhits,_maxnhits;
    int _maxarcgap;
    bool _purehits;
// vector of detector Ids corresponding to entrance and midplane
    std::vector<int> _midvids;
    std::vector<int> _entvids;
    std::vector<int> _xitvids;
// tracker ID
    int _trackerid;
// trk tuple variables
    public:
    TTree *_trkdiag;
    Int_t _fitstatus;
    Float_t _t00;
    Float_t _t00err;
    Float_t _t0;
    Float_t _t0err;
    Int_t _mcpdgid, _mcgenid, _mcproc;
    MCTrkInfo _mcinfo;
    MCTrkInfo _mcentinfo;
    MCTrkInfo _mcmidinfo;
    MCTrkInfo _mcxitinfo;
    Int_t _nhits;
    Int_t _ndof;
    Int_t _niter;
    Int_t _nt0iter;
    Int_t _nweediter;
    Int_t _nactive;
    Int_t _ncactive;
    Int_t _narcs;
    Int_t _nchits;
    Int_t _ncgood;
    Float_t _chisq;
    Float_t _fitcon;
    Float_t _radlen;
    Float_t _fitmom;
    Float_t _fitmomerr;
    Float_t _firstflt, _lastflt;
    Int_t _nsites;
    Float_t _seedmom;
    helixpar _fitpar;
    helixpar _fiterr;
    Float_t _bremsesum;
    Float_t _bremsemax;
    Float_t _bremsz;
 
    std::vector<TrkStrawHitInfo> _tshinfo;
    std::vector<TrkArcInfo> _tainfo;

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
    threevec _mchpos;
    Float_t _mcdmid;
    Float_t _mchitt0;
    Float_t _mcrdrift;
    Float_t _pdist;
    Float_t _pperp;
    Float_t _pmom;
  };
}

#endif

