//
// MC functions associated with KalFit
// $Id: KalDiag.hh,v 1.3 2014/09/20 14:34:22 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/09/20 14:34:22 $
//
#ifndef KalDiag_HH
#define KalDiag_HH
// structs
#include "TrkDiag/inc/MCEvtData.hh"
#include "TrkDiag/inc/helixpar.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawMatInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
// data
#include "art/Framework/Principal/fwd.h"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
// Utilities
#include "Mu2eUtilities/inc/MVATools.hh"
// MC data
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// MC info
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
// Mu2e tracking
#include "TrkReco/inc/TrkStrawHit.hh"
#include "TrkReco/inc/TrkDef.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "Rtypes.h"
#include "TTree.h"
#include "TClass.h"
#include "TString.h"
// C++
#include <vector>
#include <string>
#include <memory>

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
// create Branches.  The user can provide their own tree, or one will be created automatically if necessary
    TTree* createTrkDiag(TTree* trkdiag=0,const char* branchprefix="");
// find associated sim particles to a track.  The first returns a hit-weighted vector of
// all particles, the second just the one with the most hits
    void findMCTrk(const KalRep* krep,std::vector<spcount>& sct);
    void findMCTrk(const KalRep* krep,art::Ptr<SimParticle>& spp);
// access to MC data
    MCEvtData const& mcData() const { return _mcdata; }
    void findMCSteps(StepPointMCCollection const* mcsteps, cet::map_vector_key const& trkid, std::vector<int> const& vids,
	std::vector<MCStepItr>& steps);
// access to event-specific MC truth for conversion electron
    std::vector<int> const& VDids(TRACKERPOS tpos) const;
// functions to fill track information from KalRep
    void fillTrkInfo(const KalRep* krep,TrkInfo& trkinfo);
    void fillTrkFitInfo(const KalRep* krep,double fltlen,TrkFitInfo& trkfitinfo);
// MC info about a track
    void fillTrkInfoMC(art::Ptr<SimParticle> const& spp,const KalRep* krep,TrkInfoMC& trkinfomc);
    void fillTrkInfoMCStep(MCStepItr const& imcs, TrkInfoMCStep& trkinfomcstep) const;
    void fillTrkInfoMCStep(art::Ptr<SimParticle> const& spp, TrkInfoMCStep& trkinfomcstep) const; 
// hit information
    void fillHitInfo(const KalRep* krep, std::vector<TrkStrawHitInfo>& hitinfos) const;
    void fillHitInfo(const TrkStrawHit* tsh,TrkStrawHitInfo& tshinfo) const;
    void fillMatInfo(const KalRep* krep, std::vector<TrkStrawMatInfo>& hitinfos) const;
    bool fillMatInfo(const KalMaterial* ,TrkStrawMatInfo& tshinfo) const;
    void fillHitInfoMC(art::Ptr<SimParticle> const& primary,const KalRep* krep,std::vector<TrkStrawHitInfoMC>& tshinfomc) const;
    void fillHitInfoMC(art::Ptr<SimParticle> const& pspp, StrawDigiMC const& mcdigi,Straw const& straw, 
    TrkStrawHitInfoMC& tshinfomc) const;
// relationship information
    static relation relationship(art::Ptr<SimParticle> const& sppi,art::Ptr<SimParticle> const& sppj);
    static relation relationship(StrawDigiMC const& mcd1, StrawDigiMC const& mcd2);
// count CE hits
    unsigned countCEHits() const;
  private:
// cache of event data
    MCEvtData _mcdata;
// branch prefix
    std::string _branchprefix;
// event data labels
    std::string _mcptrlabel;
    std::string _mcstepslabel;
    std::string _simpartslabel, _simpartsinstance;
    std::string _mcdigislabel;
// time offsets
    SimParticleTimeOffset _toff;
// helper functions
    static void findRelatives(PtrStepPointMCVector const& mcptr,std::map<SPPtr,SPPtr>& mdmap );
    void fillTrkInfoMCStep(CLHEP::Hep3Vector const& mom, CLHEP::Hep3Vector const& pos, double charge, TrkInfoMCStep& einfo) const;
    void countHits(const KalRep* krep, TrkInfo& tinfo) const;
    void fillTrkQual(TrkInfo& trkinfo) const;
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
// track quality computation
    std::unique_ptr<MVATools> _trkqualmva;
    std::string _trkqualweights;
// struct for track info
    TrkInfo _trkinfo;
// hit information
    std::vector<TrkStrawHitInfo> _tshinfo;
    std::vector<TrkStrawMatInfo> _tminfo;
// MC true tuple variables
    TrkInfoMC _mcinfo;
    TrkInfoMCStep _mcgeninfo;
    TrkInfoMCStep _mcentinfo;
    TrkInfoMCStep _mcmidinfo;
    TrkInfoMCStep _mcxitinfo;
    std::vector<TrkStrawHitInfoMC> _tshinfomc;

  };
}

#endif

