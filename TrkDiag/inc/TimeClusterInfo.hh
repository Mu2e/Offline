//
// information about time clusters for diagnostics
// original author Dave Brown (LBNL) 8/9/2016
//
#ifndef TrkDiag_TimeClusterInfo_HH
#define TrkDiag_TimeClusterInfo_HH
// mu2e includes
#include "DataProducts/inc/XYZVec.hh"
// root includes
#include "Rtypes.h"
// C++ includes
#include <functional>
#include <string>
namespace mu2e {
  struct TimeClusterInfo {
    TimeClusterInfo() { reset(); }
    Int_t _tcindex; // cluster ID
    Int_t _nhits; // # of hits in this cluster
    Int_t _ncehits; // # of CE hits in this cluster
    Float_t _time; // cluster time
    Float_t _terr; // cluster time error
    Float_t _minhtime, _maxhtime; // min and max cluster hit time
    Float_t _maxover; // max overlap with another time cluster this event
    XYZVec _pos; // average position of cluster
    Float_t _ecalo; // calo cluster energy
    Float_t _tcalo; // calo cluster time
    Float_t _dtcalo; // calo cluster time
    XYZVec _cog; // calo cluster position
    Int_t _pripdg, _prigen, _priproc;
    Float_t _prifrac; // MC truth about the primary of a cluster
    
    void reset() { _tcindex = _prigen = _priproc = -1; _nhits = _ncehits = 0; _time = _terr = _maxover = _ecalo = _tcalo = _dtcalo = _prifrac = 0.0; _pos = _cog= XYZVec(); }
    static std::string leafnames() { 
      static std::string leaves; leaves =
      std::string("tcindex/I:nhits/I:ncehits/I:time/F:terr/F:minhtime/F:maxhtime/F:maxover/F:posx/F:posy/F:posz/F")
      + std::string(":ecalo/F:tcalo/F:dtcalo/F")
      + std::string(":cogx/F:cogy/F:cogz/F")
      + std::string(":pripdg/I:prigen/I:priproc/I:prifrac/F");
      return leaves;
    } 
  };

  struct TimeClusterHitInfo {
    TimeClusterHitInfo() { reset(); }
    Float_t _time; // hit time
    Float_t _dt; // time relative to the cluster
    Float_t _wdist; // distance along wire
    Float_t _werr; // distance error along wire
    Float_t _dphi; // resolved phi relative to the cluster
    Float_t _rho; // transverse radius of this hit
    Float_t _z; // z of this hit
    Float_t _edep; // edep
    Float_t _mva; // cluster MVA output
    Int_t _nsh, _plane;
    Int_t _mcpdg, _mcgen, _mcproc, _mcrel; // MC truth info for this hit
    Float_t _mctime, _mcmom;
    void reset() { _time = _dt = _wdist = _werr = _dphi = _rho = _z = _edep = _mva = _mctime = _mcmom = -1000.0; _mcpdg = _mcgen = _mcproc = _mcrel = 0; _nsh = _plane = -1; }
  };
    
  struct MCClusterInfo {  
    MCClusterInfo() { reset(); }
    Int_t	_nce; // # of conversion electron hits
    Int_t	_ncesel; // # of selected conversion electron hits
    Int_t	_nceclust; // # of conversion electron hits found in clusters
    Float_t	_time; // average time of CE hits (doesn't include drift!)
    XYZVec	_pos; // average position of cluste
    Float_t	_maxdphi; // max dphi WRT average
    Float_t	_minrho; // min rho WRT average
    Float_t	_maxrho; // max rho WRT average
    void reset() { _nce = _ncesel = _nceclust = 0; _time = _maxdphi = _maxrho = 0.0; _minrho = 1000.0; _pos = XYZVec();}
    static std::string leafnames() {
      static std::string leaves; leaves =
	std::string("nce/I:ncesel/I:nceclust/I:time/F")
	+std::string(":posx/F:posy/F:posz/F")
	+std::string(":maxdphi/F:minrho/F:maxrho/F");
      return leaves;
    } 
  };
// predicate to sort by decreasing # of CE hits
  struct NCEComp : public std::binary_function<TimeClusterInfo,TimeClusterInfo, bool> {
    bool operator()(TimeClusterInfo const& x,TimeClusterInfo const& y) { return x._ncehits > y._ncehits; }
  };

}
#endif
