#ifndef TrackerMC_StrawCluster_hh
#define TrackerMC_StrawCluster_hh
//
// StrawCluster represents charge arriving at one end of the straw from a single ionization.
// It does not include time folding
// or electronics effects.  It does include charge collection effects,
// such as avalanche fluctuations, trapping, gas quenching, and propagation
// effects such as drift time, wire propagation time delay, and dispersion
//
// $Id: StrawCluster.hh,v 1.4 2014/03/25 22:14:39 brownd Exp $
// $Author: brownd $
// $Date: 2014/03/25 22:14:39 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>

// Mu2e includes
#include "DataProducts/inc/StrawIndex.hh"
#include "RecoDataProducts/inc/StrawEnd.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
// toolkit includes
#include "art/Persistency/Common/Ptr.h"
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
	explicit StrawCluster(const StrawCluster& primary, StrawIndex const& index, double xfactor);
	// ghost constructor
	explicit StrawCluster(const StrawCluster& primary, double deltat);
	StrawCluster(ClusterType type,
	    StrawIndex sindex,
	    StrawEnd end,
	    double time,
	    double charge,
	    double ddist,
	    double wiredist,
	    art::Ptr<StepPointMC> const& stepmc,
	    CLHEP::HepLorentzVector const& cpos);

	StrawCluster& operator = (StrawCluster const& other);

	// Accessors
	ClusterType type() const { return _type; }
	StrawIndex strawIndex() const { return _strawIndex; }
	StrawEnd strawEnd() const { return _end; }
	double   time()       const { return _time;}
	double   charge()  const { return _charge; }
	double   driftDistance() const { return _ddist; }
	double   wireDistance() const { return _wdist; } 
	art::Ptr<StepPointMC> const&  stepPointMC() const { return _stepMC; }
	CLHEP::HepLorentzVector const& clusterPosition() const { return _cpos; }
	// Print contents of the object.
	void print( std::ostream& ost = std::cout, bool doEndl = true ) const;
      private:
	ClusterType _type; // type of clust
	StrawIndex  _strawIndex;      // Straw index
	StrawEnd	_end;		  // which end of the straw
	double     _time;            // microbunch time at the wire end, in ns
	double     _charge;          // charge at the wire end, in units of pC
	double	_ddist;		  // drift distance charge traveled to the wire
	double	_wdist;		  // distance along the wire the charge has traveled, used to calculate dispersion
	art::Ptr<StepPointMC>  _stepMC;	  // Ptr into StepPointMC collection
	CLHEP::HepLorentzVector _cpos; // position and time of the cluster that created this clust
    };
  } // namespace TrackerMC
} // namespace mu2e
#endif

