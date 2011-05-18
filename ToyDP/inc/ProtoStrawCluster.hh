#ifndef ToyDP_ProtoStrawCluster_hh
#define ToyDP_ProtoStrawCluster_hh
//
// A crude step along the way to forming real clusters.
//
// $Id: ProtoStrawCluster.hh,v 1.6 2011/05/18 15:06:33 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 15:06:33 $
//
//
#include <vector>

#include "TrackerGeom/inc/SectorId.hh"

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

#endif /* ToyDP_ProtoStrawCluster_hh */
