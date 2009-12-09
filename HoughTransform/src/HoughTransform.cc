//
// HoughTransform for circles in the L-tracker
// 
//
// $Id: HoughTransform.cc,v 1.1 2009/12/09 17:34:28 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2009/12/09 17:34:28 $
//
// Original author R.Bernstein
//
#include "HoughTransform/inc/HoughTransform.hh"
using namespace std;
namespace mu2e{
  namespace houghtransform{

  
    //    void HoughTransform::foundHoughTracks(GeomHandle<LTracker>& ltracker,edm::Handle<StepPointMCCollection>& hits,
    void HoughTransform::foundHoughTracks(GeomHandle<LTracker>& ltracker,StepPointMCCollection const* hits,
					  houghCandidates& houghCircles )
    {

      //we're looking to fit three points to a circle, so I'm going to iterate over all combinations of three hits
      StepPointMCCollection::const_iterator endhits = hits->end();

      //however, don't waste time drawing circles through the neighbor straws
      // that get passing through a module; so first compactify.  Use the HitCluster class.
      //we're going to make little clusters and then merge nearest neighbors.  
      //Any straw with no nearest neighbors will be its own cluster.

      //make the very first cluster a seed cluster for future searches.
      vector<mu2e::hitcluster::HitCluster> finalClusters;

      //and make a trial cluster based on this first hit. 
      const StepPointMC& hit = *(hits->begin());
      Straw const& straw = ltracker->getStraw( StrawIndex(hit.volumeId()) );
 
     //and now stick it into finalClusters
      mu2e::hitcluster::HitCluster initialCluster(hit,straw,hits);
      finalClusters.push_back(initialCluster);

      //loop over hits and form clusters of adjacent hits
      for (StepPointMCCollection::const_iterator ithhit = hits->begin()+1 ; ithhit < endhits; ++ithhit)
	{
	  const StepPointMC& hit = *ithhit;
	  Straw const& straw = ltracker->getStraw( StrawIndex(hit.volumeId()) );

	  //make a trial cluster based on this hit
	  mu2e::hitcluster::HitCluster trialCluster(hit,straw,hits);

	  //now loop over existing clusters. If any straws match, merge the clusters and add all straws to 
	  //the first one. go through existing clusters and compare.  I could have avoided the check here by waiting
	  //to put the first Cluster into final Clusters afterwards and keeping finalClusters.size=0,
	  //but this seems easier to follow. 

	  bool foundMatchingStraw = false;
	  for(vector<mu2e::hitcluster::HitCluster>::size_type ithCluster = 1;
	      (ithCluster <=  finalClusters.size()) && !foundMatchingStraw;
	      ++ithCluster)
	    {
	      trialCluster.matchAndMerge(foundMatchingStraw,finalClusters);
	    }
	  ///if nobody matched, this is the first time we've seen these straws and put them into
	  //the finalCluster vector.


	  if (!foundMatchingStraw) finalClusters.push_back(trialCluster);

	}

      for(vector<mu2e::hitcluster::HitCluster>::size_type ithCluster = 0; ithCluster < finalClusters.size();
	  ++ithCluster)
	{
	  //have looked through all hits and formed clusters. 
	  //the final cluster vectors may have lots of duplicates from merges and we should clean that up.
	  //it might be more elegant to do that as the clusters are formed rather than doing this here,
	  //but this way we only do it once.  
	  finalClusters.at(ithCluster).cleanUpDuplicates();

	  //do something useful -- compute the center of a cluster to hand off to the transform.
	  //getStraws produces a vector of candidates; each candidate has a set of straws and pointers to hits
	  //associated with the straws.
	  vector<mu2e::hitcluster::Candidate> candList = finalClusters[ithCluster].getStraws();
	  clusterCenters.push_back(computeClusterXYZ(candList));
	  clusterSize.push_back(candList.size());
	}

      //and we have lots of happy cluster centers.  If there are fewer than three, we can't make a circle. 

      //we'll use these to look at the vector of three-vector clusterCenters for readability
      int _xElt = 0;
      int _yElt = 1;
      int _zElt = 2;


      if (clusterCenters.size() > 3)
	{
	  for (clusterCenterVector::size_type iCluster1 = 0; iCluster1 < clusterCenters.size()-2 ; ++iCluster1)
	    {
	      double x1 = clusterCenters[iCluster1][_xElt];
	      double y1 = clusterCenters[iCluster1][_yElt];
	      for (clusterCenterVector::size_type iCluster2 = iCluster1 + 1; iCluster2 < clusterCenters.size()-2 ; ++iCluster2)
		{
		  double x2 = clusterCenters[iCluster2][_xElt];
		  double y2 = clusterCenters[iCluster2][_yElt];
		  for (clusterCenterVector::size_type iCluster3 = iCluster2 + 1; iCluster3 < clusterCenters.size() ; ++iCluster3)
		    {
		      double x3 = clusterCenters[iCluster3][_xElt];
		      double y3 = clusterCenters[iCluster3][_yElt];
		      //ok, let's make a new houghTrack but demand each cluster have two straws
		      if (clusterSize[iCluster1] >=2 && clusterSize[iCluster2] >=2 && clusterSize[iCluster3] >= 2)
			{
			  ++numberOfHoughTracks;  
			  int numberOfStrawsInThisCircle = clusterSize[iCluster1] + clusterSize[iCluster2] + clusterSize[iCluster3];
			  double radius, x0, y0, dca;
			  solveForCircle(x1,y1,x2,y2,x3,y3,radius,x0,y0,dca);
			  houghCircleStruct nextHoughCircle(radius,x0,y0,dca,numberOfStrawsInThisCircle);
			  houghCircles.push_back(nextHoughCircle);
			}
		    }
		}
	    }
	  }
	return; 
    }

    Hep3Vector HoughTransform::computeClusterXYZ(vector<mu2e::hitcluster::Candidate>& candClust)
    {
      double _averageX = 0.;
      double _averageY = 0.;
      double _clusterSize = candClust.size();

      if (candClust.empty())
	{ throw cms::Exception("RANGE") << "zero size cluster in computeClusterXYZ" ;}


      for (vector<int>::size_type ifoo = 0; ifoo < candClust.size(); ++ifoo)
	{
	  //	 	  cout << candClust.at(ifoo).id << *(candClust.at(ifoo).hitPointer) << " " ;
	  const Hep3Vector& pos = candClust.at(ifoo).hitPointer->position();
	  _averageX += pos[0];
	  _averageY += pos[1];
	}

	_averageX /= _clusterSize;
	_averageY /= _clusterSize;
      
      //we don't know the Z in this algorithm since the points are collapsed on the plane
      return Hep3Vector(_averageX,_averageY,0.);
    }


    void HoughTransform::solveForCircle(double& x1,double& y1,double& x2,double& y2,double& x3,double& y3,
								     double& radius,double& x0,double& y0, double& dca)
    //given three points, find the radius and center of the circle
    {
      radius =    (pow<2>(x1-x2) + pow<2>(y1-y2)) * (pow<2>(x1-x3) + pow<2>(y1-y3)) * (pow<2>(x2-x3) + pow<2>(y2-y3));
      radius /=     pow<2> (x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3);
      //      cout << "x1, x2,x3,y1,y2,y3 "<< x1 <<  " " << x2 << " " << x3 << " " << y1 << " " << y2 << " " << y3 << endl;
      // cout << "radius= " << sqrt(radius) << endl;
      radius = 0.5*sqrtOrThrow(radius,1.0E-06);

      //now the x- and y-center; first common term
      double denom = 2*(x3*(y1-y2) + x1*(y2-y3) + x2*(y3-y1));

      x0 = x3*x3*(y1-y2) + (x1*x1 + (y1-y2)*(y1-y3))*(y2-y3) + x2*x2*(y3-y1);
      y0 = -x2*x2*x3 + x1*x1*(x3-x2) + x3*(y1-y2)*(y1+y2) + x1*(x2*x2 - x3*x3 + y2*y2 - y3*y3) + x2*(x3*x3 - y1*y1 + y3*y3);

      x0 /= denom;
      y0 /= denom;

      //what's the dca of the circle to the origin?  
      dca = abs(sqrt(x0*x0 + y0*y0)- radius);
      //cout << "dca = "  << dca << endl;
      return;
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

    }
    return count;
  }

  }     // namespace HoughTransform
}       // namespace mu2e
