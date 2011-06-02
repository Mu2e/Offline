//
// StrawCluster 
//
// $Id: StrawCluster.cc,v 1.2 2011/06/02 22:50:54 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/02 22:50:54 $
//

// C++ includes
#include <ostream>

// Framework includes.
#include "art/Persistency/Provenance/ProductID.h"

// Mu2e includes
#include "RecoDataProducts/inc/StrawCluster.hh"

using namespace std;

namespace mu2e {

  StrawCluster::StrawCluster(std::vector<DPIndex>&  hitIndices)
  {   // Geometry info for the TTracker.
    // Get a reference to one of the T trackers.
    // Throw exception if not successful.
    _StrawHitIndices=hitIndices;
  }

  CLHEP::Hep3Vector StrawCluster::X(art::Event const & event) const
  {  
    const Tracker& tracker = getTrackerOrThrow();    
    CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
    for (size_t index =0;index<_StrawHitIndices.size();++index)
      {
	DPIndex const& junkie = _StrawHitIndices[index];
	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
	Straw str = tracker.getStraw(strawhit.strawIndex());
	const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
	pvec = pvec + mpvec;
      }
    return pvec;
  }

  double StrawCluster::Halflength(art::Event const & event) const
  {  
    double hlen = 0.0;
    const Tracker& tracker = getTrackerOrThrow();    
    for (size_t index =0;index<_StrawHitIndices.size();++index)
      {
	DPIndex const& junkie = _StrawHitIndices[index];
	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
	Straw str = tracker.getStraw(strawhit.strawIndex());
	if (str.getHalfLength()>hlen)
	  {
	    hlen=str.getHalfLength();
	  }
      }
    return hlen;
  }

  double StrawCluster::Energy(art::Event const & event) const
  {  
    double e =0.0;
    for (size_t index =0;index<_StrawHitIndices.size();++index)
      {
	DPIndex const& junkie = _StrawHitIndices[index];
	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
	e = e+strawhit.energyDep();
      }
    return e;
  }
  double StrawCluster::averageT(art::Event const & event) const
  {  
    double T =0.0;
    for (size_t index =0;index<_StrawHitIndices.size();++index)
      {
	DPIndex const& junkie = _StrawHitIndices[index];
	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
	T = T+strawhit.time();
      }
    T=T/double(_StrawHitIndices.size());
    return T;
  }
  double StrawCluster::averagedT(art::Event const & event) const
  {  
    double dT =0.0;
    for (size_t index =0;index<_StrawHitIndices.size();++index)
      {
	DPIndex const& junkie = _StrawHitIndices[index];
	StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
	dT = dT+strawhit.dt();
      }
    dT=dT/double(_StrawHitIndices.size());
    return dT;
  }
  DeviceId StrawCluster::did(art::Event const & event) const
  {   

    const Tracker& tracker = getTrackerOrThrow();    
    DPIndex const& junkie = _StrawHitIndices[0];
    StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    DeviceId did = sid.getDeviceId();
    return did;

  }
  SectorId StrawCluster::secid(art::Event const & event) const
  {   

    const Tracker& tracker = getTrackerOrThrow();    
    DPIndex const& junkie = _StrawHitIndices[0];
    StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    SectorId secid = sid.getSectorId();
    return secid;

  }
  CLHEP::Hep3Vector StrawCluster::dirX(art::Event const & event) const
  {  
    const Tracker& tracker = getTrackerOrThrow();    
    DPIndex const& junkie = _StrawHitIndices[0];
    StrawHit const& strawhit = *resolveDPIndex<StrawHitCollection>(event,junkie);
    Straw str = tracker.getStraw(strawhit.strawIndex());
    const CLHEP::Hep3Vector dirvec = str.getDirection();
    return dirvec;
  }
 
  LineSegmentPCA StrawCluster::linesegment(art::Event const& event) const
  {  
    CLHEP::Hep3Vector direction = dirX(event);
    CLHEP::Hep3Vector position  = X(event);
    double hlen  = Halflength(event);
    const CLHEP::Hep2Vector p0 =      
      CLHEP::Hep2Vector(position.getX()-hlen*direction.getX(),
			position.getY()-hlen*direction.getY());
    const CLHEP::Hep2Vector p1 = 
      CLHEP::Hep2Vector(position.getX()+hlen*direction.getX(),
			position.getY()+hlen*direction.getY());
    LineSegmentPCA linesegment(p0, p1);
    return linesegment;
  }
  

} // namespace mu2e
