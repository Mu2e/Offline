//
// MC functions associated with KalFit
// $Id: KalDiag.hh,v 1.3 2014/09/20 14:34:22 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/20 14:34:22 $
//
#ifndef KalDiag_HH
#define KalDiag_HH
// structs
#include "KalmanTests/inc/MCEvtData.hh"
#include "KalmanTests/inc/threevec.hh"
#include "KalmanTests/inc/helixpar.hh"
#include "KalmanTests/inc/MCTrkInfo.hh"
// data
#include "art/Framework/Principal/fwd.h"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// MC info
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// BaBar
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/TrkStrawHitInfo.hh"
#include "KalmanTests/inc/KalFit.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "Rtypes.h"
#include "TTree.h"
#include "TClass.h"
// C++
#include <vector>
#include <string>

namespace mu2e 
{  
 // some convenient typedefs    
  typedef art::Ptr<SimParticle> SPPtr;
 
   struct spcount {
    spcount() : _count(0) {}
    spcount(art::Ptr<SimParticle> const& spp) : _spp(spp), _count(1) {}
    void append(art::Ptr<SimParticle> const& sp) { if(sp == _spp)++_count; }
    bool operator ==(art::Ptr<SimParticle> const& sp) const { return _spp == sp; }
    art::Ptr<SimParticle> _spp;
    unsigned _count;
  };
  
  typedef StepPointMCCollection::const_iterator MCStepItr;
//  struct test : public binary_function<double,double,bool> {
//    bool operator()(double x, double y) { return x < y;}
//  };

//  Simple helper class to find MC information within collections
  
  class KalDiag {
  public:
    enum TRACKERPOS {trackerEnt=0,trackerMid, trackerExit};
    enum relation {none=-1,same,daughter,mother,sibling,udaughter,umother,usibling};
#ifndef __GCCXML__
    explicit KalDiag(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/
    virtual ~KalDiag();
// find MC data in the event.  This must be called each event, before the other functions
    bool findMCData(const art::Event& evt);
// diagnostic comparison of reconstructed tracks with MC truth
    void kalDiag(const KalRep* krep,bool fill=true);
// find associated sim particles to a track.  The first returns a hit-weighted vector of
// all particles, the second just the one with the most hits
    void findMCTrk(const KalRep* krep,std::vector<spcount>& sct);
    void findMCTrk(const KalRep* krep,art::Ptr<SimParticle>& spp);
// allow creating the trees
    TTree* createTrkDiag();
    TTree* createHitDiag();
// access to MC data
    MCEvtData const& mcData() const { return _mcdata; }
    void findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids,
	std::vector<MCStepItr>& steps);
// access to event-specific MC truth for conversion electron
    std::vector<int> const& VDids(TRACKERPOS tpos) const;
// MC info about a track
    void fillMCTrkInfo(MCStepItr const& imcs, MCTrkInfo& trkinfo) const;
    void fillMCTrkInfo(art::Ptr<SimParticle> const& spp, MCTrkInfo& einfo) const;
    void fillMCTrkInfo(art::Ptr<SimParticle> const& spp);
    static relation relationship(art::Ptr<SimParticle> const& sppi,art::Ptr<SimParticle> const& sppj);
    static relation relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2);
// MC track finder.  this function is deprecated
    bool trkFromMC(cet::map_vector_key const& trkid,TrkDef& mytrk);

  private:
// cache of event data
    MCEvtData _mcdata;
// event data labels
    std::string _mcptrlabel;
    std::string _mcstepslabel;
    std::string _simpartslabel, _simpartsinstance;
    std::string _mcdigislabel;
// time offsets
    SimParticleTimeOffset _toff;
// helper functions
    static void findRelatives(PtrStepPointMCVector const& mcptr,std::map<SPPtr,SPPtr>& mdmap );
    void hitsDiag(const KalRep* krep,art::Ptr<SimParticle> const& primary);
    const helixpar& MCHelix(TRACKERPOS tpos) const;
    void reset();
    // config parameters
    bool _fillmc;
    int _debug,_diag;
    bool _uresid;
    double _mingood;
    double _mintrkmom; // minimum true momentum at z=0 to create a track from
    double _mct0err;
    bool _mcambig;
    unsigned _minnhits,_maxnhits;
    bool _purehits;
// vector of detector Ids corresponding to entrance and midplane
    std::vector<int> _midvids;
    std::vector<int> _entvids;
    std::vector<int> _xitvids;
// trk tuple variables
    public:
    TTree *_trkdiag;
    Int_t _fitstatus;
    Float_t _t0;
    Float_t _t0err;
    Int_t _nhits;
    Int_t _ndof;
    Int_t _niter;
    Int_t _nactive;
    Int_t _ndouble,_ndactive;
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
    std::vector<TrkStrawHitInfo> _tshinfo;

// MC true tuple variables
    Int_t _npdigi, _npdgood;
    Int_t _nmc;
    Int_t _nmcactive;
    Int_t _nmchits, _nmcgood, _nmcambig;
    Int_t _mcpdgid, _mcgenid, _mcproc;
    Int_t _mcppdgid, _mcpgenid, _mcpproc;
    MCTrkInfo _mcinfo;
    MCTrkInfo _mcentinfo;
    MCTrkInfo _mcmidinfo;
    MCTrkInfo _mcxitinfo;
    std::vector<TrkStrawHitInfoMC> _tshinfomc;

  };
}

#endif

