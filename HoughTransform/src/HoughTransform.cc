//
// code for finding HoughTransform for circles in the L-tracker
// 
//
// $Id: HoughTransform.cc,v 1.5 2010/05/17 21:47:33 genser Exp $
// $Author: genser $ 
// $Date: 2010/05/17 21:47:33 $
//
// Original author R.Bernstein
//
#include "HoughTransform/inc/HoughTransform.hh"
using namespace std;
namespace mu2e{
  namespace houghtransform{

    void HoughTransform::foundHoughTracks(double radius, 
          GeomHandle<LTracker>& ltracker, houghCandidates& houghCircles )
    {
      
      //we're looking to fit all pairs of points to a circle using the 
      // known (or assumed) radius

      //we'll use these to look at the vector of three-vector clusterCenters for readability

     // need to explicitly bail if not enough clusters, to avoid 
     // -1>0 using unsigned ints...

      if (_clusterCenters.size()<2) return;

      for (clusterCenterVector::size_type iCluster1 = 0; iCluster1 < _clusterCenters.size()-2 ; ++iCluster1)
      {
        for (clusterCenterVector::size_type iCluster2 = iCluster1 + 1; iCluster2 < _clusterCenters.size()-1 ; ++iCluster2)
        {
             //ok, let's make a new houghTrack but demand each cluster have two straws
             if (_clusterSizes[iCluster1] >=2 && _clusterSizes[iCluster2] >=2 )
             {
               ++numberOfHoughTracks;  
               int numberOfStrawsInThisCircle = _clusterSizes[iCluster1] + 
                                                _clusterSizes[iCluster2];
               double x0m, y0m, dcam, x0p, y0p, dcap;
               if (solveForCircle2P(_clusterCenters[iCluster1],
                                _clusterCenters[iCluster2],
                                radius,x0m,y0m,dcam,x0p,y0p,dcap)) {
                  // we get two solutions when using 2 points+radius
                  houghCircleStruct nextHoughCirclem(radius,x0m,y0m,dcam,
                                                   numberOfStrawsInThisCircle);
                  houghCircles.push_back(nextHoughCirclem);
                  houghCircleStruct nextHoughCirclep(radius,x0p,y0p,dcap,
                                                   numberOfStrawsInThisCircle);
                  houghCircles.push_back(nextHoughCirclep);
               }// solution found

             } // enough straws

        } // cluster 2

      } // cluster 1

      return; 

    } //foundHoughTracks(2 points + radius,...)
  
    void HoughTransform::foundHoughTracks(GeomHandle<LTracker>& ltracker,
                  houghCandidates& houghCircles )
    {

      //we're looking to fit three points to a circle, so I'm going to iterate over all combinations of three hits

      //however, don't waste time drawing circles through the neighbor straws
      // that get passing through a module; so first compactify.  Use the HitCluster class.
      //we're going to make little clusters and then merge nearest neighbors.  
      //Any straw with no nearest neighbors will be its own cluster.

/*

    Deleted a bunch of stuff about making HitClusters, and moved it to 
    its own function to be called by the plugin, since eventually, HitClusters 
    will be created by a module upstream of the Hough finding algorithm.

*/

      //we'll use these to look at the vector of three-vector clusterCenters for readability

     // need to explicitly bail if not enough clusters, to avoid 
     // -1>0 using unsigned ints...

      if (_clusterCenters.size()<3) return;

      for (clusterCenterVector::size_type iCluster1 = 0; iCluster1 < _clusterCenters.size()-2 ; ++iCluster1)
      {
        for (clusterCenterVector::size_type iCluster2 = iCluster1 + 1; iCluster2 < _clusterCenters.size()-1 ; ++iCluster2)
        {
          for (clusterCenterVector::size_type iCluster3 = iCluster2 + 1; iCluster3 < _clusterCenters.size() ; ++iCluster3)
          {
             //ok, let's make a new houghTrack but demand each cluster have two straws
             if (_clusterSizes[iCluster1] >=2 && 
             _clusterSizes[iCluster2] >=2 && 
             _clusterSizes[iCluster3] >= 2)
             {
               ++numberOfHoughTracks;  
               int numberOfStrawsInThisCircle = _clusterSizes[iCluster1] + 
                                                _clusterSizes[iCluster2] + 
                                                _clusterSizes[iCluster3];
               double radius, x0, y0, dca;
               bool found=false;
               if (_use3P) {
                 //geometric version
                 found=solveForCircle3P(_clusterCenters[iCluster1],
                                        _clusterCenters[iCluster2],
                                        _clusterCenters[iCluster3], 
                                                        radius,x0,y0,dca);
               } else {
                 //algrebraic version
                 found=solveForCircle(_clusterCenters[iCluster1],
                                      _clusterCenters[iCluster2],
                                      _clusterCenters[iCluster3], 
                                                        radius,x0,y0,dca);
               }// which type of circle solver

               if (found) {
                  houghCircleStruct nextHoughCircle(radius,x0,y0,dca,
                                                numberOfStrawsInThisCircle);
                  houghCircles.push_back(nextHoughCircle);
               } //solution found

             } // enough straws

          } // cluster 3

        } // cluster 2

      } // cluster 1

      return; 

    } //foundHoughTracks(3 points)

    CLHEP::Hep3Vector HoughTransform::computeClusterXYZ(vector<mu2e::hitcluster::Candidate>& candClust)
    {
      double averageX = 0.;
      double averageY = 0.;
      double clusterSize = candClust.size();

      if (candClust.empty())
	{ throw cms::Exception("RANGE") << "zero size cluster in computeClusterXYZ" ;}


      for (vector<int>::size_type ifoo = 0; ifoo < candClust.size(); ++ifoo)
	{
	  //	 	  cout << candClust.at(ifoo).id << *(candClust.at(ifoo).hitPointer) << " " ;
	  const CLHEP::Hep3Vector& pos = candClust.at(ifoo).hitPointer->position();
	  averageX += pos[0];
	  averageY += pos[1];
	}

	averageX /= clusterSize;
	averageY /= clusterSize;
      
      //we don't know the Z in this algorithm since the points are collapsed on the plane
      return CLHEP::Hep3Vector(averageX,averageY,0.);
    }


    bool HoughTransform::solveForCircle2P(const CLHEP::Hep3Vector& v1,
                                          const CLHEP::Hep3Vector& v2,
                     double radius,
                     double& x0m, double& y0m, double& dcam,
                     double& x0p, double& y0p, double& dcap)
    //given two points and the radius, find one of the two solutions for 
    //the center. returns false for no solution.  sign<0 for - , all else+
    {

      CLHEP::Hep2Vector p1(v1.x(),v1.y());
      CLHEP::Hep2Vector p2(v2.x(),v2.y());

      // vector between the two points (the chord)
      CLHEP::Hep2Vector p1p2(p2-p1);

      // distance p1p2
      double h=p1p2.mag();     

      // check: are points close enough to be on circle of given radius
      if (h>radius*2.) return false;

      // make it a unit vector
      p1p2.setMag(1);

      // unit vector perpendicular to chord
      CLHEP::Hep2Vector perp(-p1p2.y(),p1p2.x());

      // midpoint of the chord
      CLHEP::Hep2Vector m(p1+p2); m*=0.5;

      // distance from midpoint of the chord to the center of the circle
      double d=sqrt(radius*radius-h*h/4.);

      CLHEP::Hep2Vector center;
      CLHEP::Hep2Vector centerm=m-perp*d; // "-" solution
      CLHEP::Hep2Vector centerp=m+perp*d; // "+" solution

      // set the return values
      x0m=centerm.x();
      y0m=centerm.y();
      dcam=TMath::Abs(centerm.mag()-radius);
      x0p=centerp.x();
      y0p=centerp.y();
      dcap=TMath::Abs(centerp.mag()-radius);

      return true;
    }

    bool HoughTransform::solveForCircle3P(const CLHEP::Hep3Vector& v1,
                                          const CLHEP::Hep3Vector& v2,
                                          const CLHEP::Hep3Vector& v3,
                          double& radius,double& x0,double& y0, double& dca)
    //given three points, find the radius and center of the circle using 
    // geometric approach
    {

      // perpendicular bisector of any two points is a radial
      
      // intersection of any two of those is the center

      // this assumes we only care about x-y plane, but we use 3 vectors
      // to get the cross product.  Force z=0 on input points.

      CLHEP::Hep3Vector p1(v1); p1.setZ(0);
      CLHEP::Hep3Vector p2(v2); p2.setZ(0);
      CLHEP::Hep3Vector p3(v3); p3.setZ(0);

      CLHEP::Hep3Vector uz(0,0,1);// unit z vector

      // 2 vs. 1
      CLHEP::Hep3Vector c12(p1-p2); // chord
      CLHEP::Hep3Vector m12(p1+p2); m12/=2; // midpoint of the chord
      CLHEP::Hep3Vector p12=c12.cross(uz); // parallel to the radial through
                                           // the chord

      CLHEP::Hep3Vector c23(p2-p3); // chord
      CLHEP::Hep3Vector m23(p2+p3); m23/=2; // midpoint of the chord
      CLHEP::Hep3Vector p23=c23.cross(uz); // parallel to the radial through
                                           // the chord

      //center is at intersection of two lines, which are parameterized by
      //  m12+s*p12 and m23+t*p23.  s and t are scalars
      // set them equal, and take the cross product with p12.  That kills s.
      // then take the dot product with m12xp12 to get a scalar expression
      // for t: 
      // (m12xp12).(m12xp12) = (m23xp12).(m12xp12) + t * (p23xp12).(m12xp12).  
      // Plug back in to find center, then radius, dca
      CLHEP::Hep3Vector m12p12=m12.cross(p12);
      CLHEP::Hep3Vector m23p12=m23.cross(p12);
      CLHEP::Hep3Vector p23p12=p23.cross(p12);

      double denom=p23p12.dot(m12p12);
      if (TMath::Abs(denom)<1e-10) return false;

      double t=(m12p12.dot(m12p12)-m23p12.dot(m12p12))/denom;
      CLHEP::Hep3Vector center(p23); center*=t; center+=m23;
      x0=center.x(); y0=center.y();

      CLHEP::Hep3Vector radial=center; radial-=p1;
      radius=radial.mag();
      
      dca=TMath::Abs(center.mag()-radius);
      
      return true;
    }

    bool HoughTransform::solveForCircle(const CLHEP::Hep3Vector& v1,
                                          const CLHEP::Hep3Vector& v2,
                                          const CLHEP::Hep3Vector& v3,
                          double& radius,double& x0,double& y0, double& dca)
    //given three points, find the radius and center of the circle using
    // algebraic approach
    {

      double x1=v1.x();
      double y1=v1.y();
      double x2=v2.x();
      double y2=v2.y();
      double x3=v3.x();
      double y3=v3.y();

      radius =    (pow<2>(x1-x2) + pow<2>(y1-y2)) * 
                  (pow<2>(x1-x3) + pow<2>(y1-y3)) * 
                  (pow<2>(x2-x3) + pow<2>(y2-y3));
      radius /=     pow<2> (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3);
      radius = 0.5*sqrtOrThrow(radius,1.0E-06);

      //now the x- and y-center; first common term
      double denom = 2*(x3*(y1-y2) + x1*(y2-y3) + x2*(y3-y1));

      if (TMath::Abs(denom)<1e-10) return false;// protect against /0

      x0 = x3*x3*(y1-y2) + (x1*x1 + (y1-y2)*(y1-y3))*(y2-y3) + x2*x2*(y3-y1);
      y0 = -x2*x2*x3 + x1*x1*(x3-x2) + x3*(y1-y2)*(y1+y2) + x1*(x2*x2 - x3*x3 + y2*y2 - y3*y3) + x2*(x3*x3 - y1*y1 + y3*y3);

      x0 /= denom;
      y0 /= denom;

      //what's the dca of the circle to the origin?  
      dca = abs(sqrt(x0*x0 + y0*y0)- radius);

      return true;
    }

    int HoughTransform::countHitNeighbours( Straw const& straw, 
					    //				    edm::Handle<StepPointMCCollection>& hits ){
				    StepPointMCCollection const* hits ){
    
      int count(0);
      vector<StrawIndex> const& nearest = straw.nearestNeighboursByIndex();
      for ( vector<int>::size_type ihit =0;
	    ihit<nearest.size(); ++ihit ){
  
        StrawIndex idx = nearest[ihit];
  
        for( StepPointMCCollection::const_iterator 
	       i = hits->begin(),
	       e = hits->end(); i!=e ; ++i ) {
	  const StepPointMC& hit = *i;
	  if ( hit.volumeId() == idx.asInt() ){
	    ++count;
	    break;
	  }
        }

      }// ihit
      return count;
    }// countHitNeighbors

    void HoughTransform::FindCenters()
    {
      for(ClusterList::size_type iclus = 0; iclus < _hitClusters.size();
	  ++iclus)
	{

	  //do something useful -- compute the center of a cluster to hand off to the transform.
	  //getStraws produces a vector of candidates; each candidate has a set of straws and pointers to hits
	  //associated with the straws.
	  vector<mu2e::hitcluster::Candidate> candList = 
                              _hitClusters[iclus].getStraws();
	  _clusterCenters.push_back(computeClusterXYZ(candList));
	  _clusterSizes.push_back(candList.size());
	}

    } // FindCenters

    /*static*/ void HoughTransform::MakeClusters(StepPointMCCollection const* hits,
                           ClusterList& newClusters)
    {

      // Master geometry for the LTracker.
      GeomHandle<LTracker> ltracker;

      //and make a trial cluster based on first hit. 
      const StepPointMC& hit = *(hits->begin());
      Straw const& straw = ltracker->getStraw( StrawIndex(hit.volumeId()) );
 
     //and now stick it into newClusters
      mu2e::hitcluster::HitCluster initialCluster(hit,straw,hits);
      newClusters.push_back(initialCluster);

      StepPointMCCollection::const_iterator ithhit = hits->begin();
      const StepPointMCCollection::const_iterator endhit = hits->end();

      //loop over hits and form clusters of adjacent hits
      for ( ; ithhit < endhit; ++ithhit)
	{
	  const StepPointMC& hit = *ithhit;
	  Straw const& straw = ltracker->getStraw( StrawIndex(hit.volumeId()) );

	  //make a trial cluster based on this hit
	  mu2e::hitcluster::HitCluster trialCluster(hit,straw,hits);

	  //now loop over existing clusters. If any straws match, merge the clusters and add all straws to 
	  //the first one. go through existing clusters and compare.  I could have avoided the check here by waiting
	  //to put the first Cluster into newClusters Clusters afterwards and keeping newClusters.size=0,
	  //but this seems easier to follow. 

	  bool foundMatchingStraw = false;
	  for(vector<mu2e::hitcluster::HitCluster>::size_type ithCluster = 1;
	      (ithCluster <=  newClusters.size()) && !foundMatchingStraw;
	      ++ithCluster)
	    {
	      trialCluster.matchAndMerge(foundMatchingStraw,newClusters);
	    }
	  ///if nobody matched, this is the first time we've seen these straws and put them into
	  //the newClusters vector.


	  if (!foundMatchingStraw) newClusters.push_back(trialCluster);

	}//hit

        for(ClusterList::size_type ithCluster = 0; 
                            ithCluster < newClusters.size(); ++ithCluster)
	{
	  //have looked through all hits and formed clusters. 
	  //the newClusters cluster vectors may have lots of duplicates from merges and we should clean that up.
	  //it might be more elegant to do that as the clusters are formed rather than doing this here,
	  //but this way we only do it once.  
	  newClusters.at(ithCluster).cleanUpDuplicates();

	}

    }// MakeClusters

  }     // namespace HoughTransform
}       // namespace mu2e
