#ifndef TrackerMC_StrawClusterSequence_hh
#define TrackerMC_StrawClusterSequence_hh
//
// StrawClusterSequence is a time-ordered sequence of StrawClusters
//
// $Id: StrawClusterSequence.hh,v 1.1 2013/12/07 19:51:42 brownd Exp $
// $Author: brownd $
// $Date: 2013/12/07 19:51:42 $
//
// Original author David Brown, LBNL
//

// C++ includes
#include <iostream>
#include <list>
// Mu2e includes
#include "TrackerMC/inc/StrawCluster.hh"
#include "DataProducts/inc/StrawIndex.hh"

namespace mu2e {
  typedef std::list<StrawCluster> ClusterList;
  class StrawClusterSequence {
    public:
// constructors
      StrawClusterSequence();
      StrawClusterSequence(StrawCluster const& clust);
      StrawClusterSequence(StrawIndex const& index, StrawEnd end);
      StrawClusterSequence(StrawClusterSequence const& other);
      StrawClusterSequence& operator =(StrawClusterSequence const& other);
      // accessors: just hand over the list!
      ClusterList const& clustList() const { return _clist; }
      // insert a new clust, in time order.
      ClusterList::iterator insert(StrawCluster const& clust);
      StrawIndex strawIndex() const { return _strawIndex; }
      StrawEnd strawEnd() const { return _end; }
    private:
      StrawIndex _strawIndex;
      StrawEnd _end;
      ClusterList _clist; // time-ordered sequence of clusts
  };
}
#endif


