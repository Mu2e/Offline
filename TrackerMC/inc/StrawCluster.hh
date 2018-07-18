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
#include "MCDataProducts/inc/StepPointMC.hh"
// toolkit includes
#include "canvas/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {
  namespace TrackerMC {
    class StrawCluster{
    public:
      enum ClusterType {unknown=-1,primary=0,xtalk=1,noise=2};
      // constructors
      StrawCluster();
      StrawCluster(const StrawCluster&);
      // x-talk constructor
      explicit StrawCluster(const StrawCluster& primary, StrawId const& sid, double xfactor);
      // ghost constructor
      explicit StrawCluster(const StrawCluster& primary, double deltat);
      StrawCluster(ClusterType type,
                   StrawId sid,
                   StrawEnd end,
                   double time,
                   double charge,
                   double ddist,
                   double phi, //JB: added
                   double wiredist,
                   double drifttime,
                   double proptime,
                   art::Ptr<StepPointMC> const& stepmc,
                   CLHEP::HepLorentzVector const& cpos);

      StrawCluster& operator = (StrawCluster const& other);

      // Accessors
      ClusterType type() const { return _type; }
      StrawId strawId() const { return _strawId; }
      StrawEnd strawEnd() const { return _end; }
      double   time()       const { return _time;}
      double   charge()  const { return _charge; }
      double   driftDistance() const { return _ddist; }
      double   phi() const { return _phi; } //JB: added
      double   wireDistance() const { return _wdist; }
      double   driftTime() const { return _drifttime; }
      double   propTime() const { return _proptime; }
      art::Ptr<StepPointMC> const&  stepPointMC() const { return _stepMC; }
      CLHEP::HepLorentzVector const& clusterPosition() const { return _cpos; }
      // Print contents of the object.
      void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
    private:
      ClusterType _type; // type of clust
      StrawId  _strawId;      // Straw id
      StrawEnd	_end;		  // which end of the straw
      double     _time;            // microbunch time at the wire end, in ns since EventWindowMarker
      double     _charge;          // charge at the wire end, in units of pC
      double	_ddist;		  // drift distance charge traveled to the wire
      double _phi;    //JB: angle between E and B at ionization event
      double	_wdist;		  // distance along the wire the charge has traveled, used to calculate dispersion
      double _drifttime;
      double _proptime;
      art::Ptr<StepPointMC>  _stepMC;	  // Ptr into StepPointMC collection
      CLHEP::HepLorentzVector _cpos; // position and time of the cluster that created this clust
    };
  } // namespace TrackerMC
} // namespace mu2e
#endif
