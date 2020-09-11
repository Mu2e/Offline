//
//
// Original author G. Pezzullo
//

#include <ostream>
#include "RecoDataProducts/inc/CaloCluster.hh"


namespace mu2e {

    void CaloCluster::print( std::ostream& ost) const 
    {
      ost << "CaloCluster :   "
          << " section: "            << diskId_
          << " time: "               << time_                 << "\n"
          << " energyDep: "          << energyDep_            << "\n"
          << " COG3Vector.u: "       << cog3Vector_.x()
          << " COG3Vector.v: "       << cog3Vector_.y()
          << " COG3Vector.w: "       << cog3Vector_.z()       <<"\n"
          << " size: "               << caloCrystalHitsPtrVector_.size()<<"\n";
    }

} 
