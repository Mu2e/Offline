#ifndef TrackerMC_StrawCluster_hh
#define TrackerMC_StrawCluster_hh
//
// StrawCluster represents charge arriving at one end of the straw from a single ionization.
// It does not include time folding
// or electronics effects.  It does include charge collection effects,
// such as avalanche fluctuations, trapping, gas quenching, and propagation
// effects such as drift time, wire propagation time delay, and dispersion
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "TrackerMC/inc/StrawPosition.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
// toolkit includes
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {
  namespace TrackerMC {
    class StrawCluster{
    public:
      enum ClusterType {unknown=-1,primary=0,xtalk=1,noise=2};
      // constructors
      StrawCluster();
      // x-talk constructor
      explicit StrawCluster(const StrawCluster& primary, StrawId const& sid, float xfactor);
      // ghost constructor
      explicit StrawCluster(const StrawCluster& primary, double deltat);
      explicit StrawCluster(ClusterType type,StrawId sid,
	  StrawEnd end,
	  float time,
	  float charge,
	  float wdist,
	  StrawPosition const& pos,
	  float drifttime,
	  float proptime,
	  art::Ptr<StrawGasStep> const& sgs,
	  float ctime);
      // use compiler version of copy, assignment
      // Accessors
      ClusterType type() const { return _type; }
      StrawId strawId() const { return _strawId; }
      StrawEnd strawEnd() const { return _end; }
      double time()       const { return _time;}
      float   charge()  const { return _charge; }
      StrawPosition const& cluPos() const { return _pos; }
      float wireDistance() const { return _wdist; }
      float driftDistance() const { return _pos.Rho(); }
      float driftPhi() const { return _pos.Phi(); }
      float   driftTime() const { return _drifttime; }
      float   propTime() const { return _proptime; }
      art::Ptr<StrawGasStep> const& strawGasStep() const { return _sgsptr; }
      float cluTime() const { return _ctime; }
      // Print contents of the object.
      void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
    private:
      ClusterType _type; // type of cluster
      StrawId  _strawId;      // Straw id
      StrawEnd	_end;		  // which end of the straw
      float  _time;            // microbunch time at the wire end, in ns since EventWindowMarker, offsets and wrapping applied
      float  _charge;          // charge at the wire end, in units of pC
      float _wdist;    // propagation distance from cluster to the wire end
      StrawPosition _pos;  // cluster position WRT the straw 
      float _drifttime; // drift time to the wire
      float _proptime;  // propagation time to the wire end
      art::Ptr<StrawGasStep> _sgsptr; 
      float _ctime; 
    };
  } // namespace TrackerMC
} // namespace mu2e
#endif
