#ifndef ProtoStrawCluster_hh
#define ProtoStrawCluster_hh
//
// A crude step along the way to forming real clusters.
//
// $Id: ProtoStrawCluster.hh,v 1.1 2009/10/28 14:14:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/28 14:14:13 $
//
//
#include <vector>

#include "LTrackerGeom/inc/SectorId.hh"

namespace mu2e {

  struct ProtoStrawCluster{

    SectorId id;
    std::vector<int> hitIndices;

    ProtoStrawCluster( SectorId const& sectorId, int index);
    
    void add( int idx){
      hitIndices.push_back(idx);
    }

    int size() const{ return hitIndices.size(); }

    //int operator[]( int i) const { return hitIndices[i]; }

    int at( int i) const { return hitIndices.at(i); }
    
  };

}

#endif
