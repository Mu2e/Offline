//
// Loop over straws added in the previous iteration and add
// their nearest neighbours to the list.
// 
// If a hit straw appears more than once in input list, then
// all of those hits to the cluster.
//
// $Id: growCluster.cc,v 1.1 2009/10/28 14:14:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/28 14:14:13 $
// 

// C++ includes
#include <iostream>

// Mu2e includes
#include "HitMakers/inc/growCluster.hh"
#include "ToyDP/inc/ProtoStrawCluster.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/CrudeStrawHitCollection.hh"

using namespace std;

namespace mu2e {
  
  int growCluster ( ProtoStrawCluster&              cluster,
		    int                             startCluster,
		    int                             startHit,
		    edm::Handle<CrudeStrawHitPData> pdataHandle,
		    std::vector<int>&               used,
		    LTracker const&                 ltracker
		    ){

    // Alias for readability.
    CrudeStrawHitPData const& hits(*pdataHandle);

    // Initialize the return value.
    int nadded(0);

    // Size of the cluster on input.
    int last(cluster.size());

    // Number of hits.
    int const nhits(hits.size());

    // Loop over hits added in the previous iteration.
    for ( int i=startCluster; i<last; ++i ){

      CrudeStrawHit const& baseHit(hits.at(cluster.at(i)));
      Straw const&         baseStraw(ltracker.getStraw(baseHit.strawIndex));

      // Loop over remaining  hits in the hit list.
      for ( int j=startHit; j<nhits; ++j ){

	// Skip hits that are already used.
	if ( j       == i ) continue;
	if ( used[j] == 1 ) continue;
	
	CrudeStrawHit const& hit(hits.at(j));
	
	// Add neighbours to the cluster and mark the hit as used.
	if ( baseStraw.isNearestNeighbour(hit.strawIndex) ) {
	  used[j] = 1;
	  cluster.add(j);
	  ++nadded;
	}
	// If the straw itself appears a second time, add it too.
	else if ( baseStraw.Index() == hit.strawIndex ){
	  used[j] = 1;
	  cluster.add(j);
	  ++nadded;
	}

      }
      
    }
    return nadded;
  }
}
