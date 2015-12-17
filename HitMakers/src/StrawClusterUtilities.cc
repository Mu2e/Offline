//
// StrawClusterUtilities
//

// C++ includes
//#include <ostream>

// Mu2e includes
#include "HitMakers/inc/StrawClusterUtilities.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"

using namespace std;

namespace mu2e {


  CLHEP::Hep3Vector StrawClusterUtilities::midX(StrawCluster const & cluster,art::Event const & event) const
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

  CLHEP::Hep3Vector StrawClusterUtilities::dirX(StrawCluster const & cluster,art::Event const & event) const
  {
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];
    Straw str = tracker.getStraw(strawhit.strawIndex());
    const CLHEP::Hep3Vector dirvec = str.getDirection();
    return dirvec;
  }

 CLHEP::Hep3Vector StrawClusterUtilities::dTX(StrawCluster const & cluster,art::Event const & event) const
  {
    const Tracker& tracker = getTrackerOrThrow();
    CLHEP::Hep3Vector pvec = CLHEP::Hep3Vector(0.,0.,0.);
    double averagedt = 0.0;
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    for (size_t index =0;index<strawHits.size();++index)
      {
        StrawHit const& strawhit = *strawHits[index];
        Straw str = tracker.getStraw(strawhit.strawIndex());
        const CLHEP::Hep3Vector mpvec  = str.getMidPoint();
        averagedt = averagedt+ strawhit.dt();
        pvec = pvec + mpvec;
      }

    double a          = 1./double(strawHits.size());
    averagedt         = averagedt * a;
    pvec              = pvec*a;
    double _timetodist= 149.8962;
    double disttomid  = averagedt* _timetodist;
    cout << "disttomid"<<disttomid<<endl;
    cout << "Halflength: "<<Halflength(cluster,event)<<endl;
    cout << "pvec: "<<pvec<<endl;
    cout << "dirX: " <<dirX(cluster,event)<<endl;
    //    CLHEP::Hep3Vector dirvec = dirX(cluster,event);
    pvec              = pvec+disttomid*dirX(cluster,event);
    cout << "pvec: "<<pvec<<endl;
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
  PlaneId StrawClusterUtilities::did(StrawCluster const & cluster,art::Event const & event) const
  {
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    PlaneId did = sid.getPlaneId();
    return did;

  }
  PanelId StrawClusterUtilities::secid(StrawCluster const & cluster,art::Event const & event) const
  {
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    PanelId secid = sid.getPanelId();
    return secid;
  }
  PlaneId StrawClusterUtilities::Station(StrawCluster const & cluster,art::Event const & event) const
  {
    const Tracker& tracker = getTrackerOrThrow();
    StrawHitPtrVector const & strawHits = cluster.strawHits();
    StrawHit const& strawhit = *strawHits[0];
    Straw    str = tracker.getStraw(strawhit.strawIndex());
    StrawId  sid = str.id();
    PlaneId did = sid.getPlaneId();
    int station = (int)(did*0.5);
    return station;

  }


  LineSegmentPCA StrawClusterUtilities::linesegment(StrawCluster const & cluster,art::Event const& event) const
  {
    CLHEP::Hep3Vector direction = dirX(cluster,event);
    CLHEP::Hep3Vector position  = midX(cluster,event);
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
  multimap<int,StrawCluster> StrawClusterUtilities::clusterbydid(StrawClusterCollection const& clusters,art::Event const& event) const
  {
    multimap<int,StrawCluster> clubydid;
    for ( size_t cluster=0; cluster<clusters.size(); ++cluster) // Loop over StrawClusters
      {
        StrawCluster const& scluster = clusters.at(cluster);
        clubydid.insert(pair<int,StrawCluster>(did(scluster,event),scluster));
      }
    return clubydid;
  }
  multimap<int,StrawCluster> StrawClusterUtilities::clusterbystation(StrawClusterCollection const& clusters,art::Event const& event) const
  {
    multimap<int,StrawCluster> clubystation;
    for ( size_t cluster=0; cluster<clusters.size(); ++cluster) // Loop over StrawClusters
      {
        StrawCluster const& scluster = clusters.at(cluster);
        clubystation.insert(pair<int,StrawCluster>(Station(scluster,event),scluster));
      }
    return clubystation;
  }
} // namespace mu2e
