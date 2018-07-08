#ifndef TrackerMC_StrawClusterSequence_hh
#define TrackerMC_StrawClusterSequence_hh
//
// StrawClusterSequence is a time-ordered sequence of StrawClusters
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <list>
// Mu2e includes
#include "TrackerMC/inc/StrawCluster.hh"
#include "DataProducts/inc/StrawId.hh"

namespace mu2e {
  namespace TrackerMC {
    typedef std::list<StrawCluster> ClusterList;
    class StrawClusterSequence {
      public:
	// constructors
	StrawClusterSequence();
	StrawClusterSequence(StrawCluster const& clust);
	StrawClusterSequence(StrawId const& sid, StrawEnd end);
	StrawClusterSequence(StrawClusterSequence const& other);
	StrawClusterSequence& operator =(StrawClusterSequence const& other);
	// accessors: just hand over the list!
	ClusterList const& clustList() const { return _clist; }
	// insert a new clust, in time order.
	ClusterList::iterator insert(StrawCluster const& clust);
	StrawId const& strawId() const { return _strawId; }
	StrawEnd const& strawEnd() const { return _end; }
      private:
	StrawId _strawId;
	StrawEnd _end;
	ClusterList _clist; // time-ordered sequence of clusts
    };
  }
}
#endif
