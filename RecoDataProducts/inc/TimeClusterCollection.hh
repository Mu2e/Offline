#ifndef RecoDataProducts_TimeClusterCollection_hh
#define RecoDataProducts_TimeClusterCollection_hh

//
// Define a type for a collection of TimeCluster objects.
//
// $Id: TimeClusterCollection.hh,v 1.0 2018/10/11 11:02 S. Middleton Exp $
// $Author: Sophie Middleton
// $Date: 2018/10/11 11:02:11 $
//
// Original author Sophie Middleton
//A Time Cluster is just a time-selected collection of StrawHits....
#include <vector>
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

namespace mu2e {
  typedef std::vector<mu2e::TimeCluster> TimeClusterCollection;
}

#endif /* RecoDataProducts_TimeClusterCollection_hh */
