#ifndef ToyDP_ProtoStrawCluster_hh
#define ToyDP_ProtoStrawCluster_hh
//
// A crude step along the way to forming real clusters.
//
// $Id: ProtoStrawCluster.hh,v 1.4 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
//
#include <vector>

#include "TrackerGeom/inc/SectorId.hh"

namespace mu2e {

  struct ProtoStrawCluster{

    SectorId id;
    std::vector<int32_t> hitIndices;

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
