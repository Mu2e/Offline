//
// StrawClusterUtilities 
//
// $Id: StrawClusterUtilities.cc,v 1.1 2011/06/09 21:21:52 wenzel Exp $
// $Author: wenzel $
// $Date: 2011/06/09 21:21:52 $
//

// C++ includes
//#include <ostream>

// Framework includes.
#include "art/Persistency/Provenance/ProductID.h"

// Mu2e includes
#include "Mu2eUtilities/inc/StrawClusterUtilities.hh"

using namespace std;

namespace mu2e {


  CLHEP::Hep3Vector StrawClusterUtilities::X(StrawCluster const & cluster,art::Event const & event) const
  {  
    const Tracker& tracker = getTrackerOrThrow();    
    CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    for (size_t index =0;index<strawHits.size();++index)
      {
	StrawHit const& strawhit = *strawHits[index];
	Straw str = tracker.getStraw(strawhit.strawIndex());
	const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
	pvec = pvec + mpvec;
      }
      double a = 1./double(strawHits.size());
      pvec = pvec*a;
    return pvec;
  }

  double StrawClusterUtilities::Halflength(StrawCluster const & cluster,art::Event const & event) const
  {  
    double hlen = 0.0;
    const Tracker& tracker = getTrackerOrThrow(); 
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    for (size_t index =0;index<strawHits.size();++index)
      {
	StrawHit const& strawhit = *strawHits[index];
	Straw str = tracker.getStraw(strawhit.strawIndex());
	if (str.getHalfLength()>hlen)
	  {
	    hlen=str.getHalfLength();
	  }
      }
    return hlen;
  }

  double StrawClusterUtilities::Energy(StrawCluster const & cluster,art::Event const & event) const
  {  
    double e =0.0;
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    for (size_t index =0;index<strawHits.size();++index)
      {
	StrawHit const& strawhit = *strawHits[index];
	e = e+strawhit.energyDep();
      }
    return e;
  }
  double StrawClusterUtilities::averageT(StrawCluster const & cluster,art::Event const & event) const
  {  
    double T =0.0;
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    for (size_t index =0;index<strawHits.size();++index)
      {
	StrawHit const& strawhit = *strawHits[index];
	T = T+strawhit.time();
      }
    T=T/double(strawHits.size());
    return T;
  }
  double StrawClusterUtilities::averagedT(StrawCluster const & cluster,art::Event const & event) const
  {  
    double dT =0.0;
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    for (size_t index =0;index<strawHits.size();++index)
      {
	StrawHit const& strawhit = *strawHits[index];
	dT = dT+strawhit.dt();
      }
    dT=dT/double(strawHits.size());
    return dT;
  }
  DeviceId StrawClusterUtilities::did(StrawCluster const & cluster,art::Event const & event) const
  {       
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];    
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    DeviceId did = sid.getDeviceId();
    return did;

  }
  SectorId StrawClusterUtilities::secid(StrawCluster const & cluster,art::Event const & event) const
  {   
    const Tracker& tracker = getTrackerOrThrow();    
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0]; 
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    SectorId secid = sid.getSectorId();
    return secid;
  }
  DeviceId StrawClusterUtilities::Station(StrawCluster const & cluster,art::Event const & event) const
  {       
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];    
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    DeviceId did = sid.getDeviceId();
    int station = (int)(did*0.5);
    return station;

  }
  CLHEP::Hep3Vector StrawClusterUtilities::dirX(StrawCluster const & cluster,art::Event const & event) const
  {  
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];
    Straw str = tracker.getStraw(strawhit.strawIndex());
    const CLHEP::Hep3Vector dirvec = str.getDirection();
    return dirvec;
  }
 
  LineSegmentPCA StrawClusterUtilities::linesegment(StrawCluster const & cluster,art::Event const& event) const
  {  
    CLHEP::Hep3Vector direction = dirX(cluster,event);
    CLHEP::Hep3Vector position  = X(cluster,event);
    double hlen  = Halflength(cluster,event);
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
