#include "ToyDP/inc/ProtoStrawCluster.hh"

namespace mu2e {

  ProtoStrawCluster::ProtoStrawCluster( SectorId const& sectorId, int index):
    id(sectorId),
    hitIndices(1,index){
  }
  
}
