//
// $Id: CaloCluster.cc,v 1.3 2012/09/06 19:58:26 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/09/06 19:58:26 $
//
// Original author G. Pezzullo
//

#include <ostream>
#include "RecoDataProducts/inc/CaloCluster.hh"


namespace mu2e {

    void CaloCluster::print( std::ostream& ost) const 
    {
      ost << "CaloCluster :   "
          << " section: "            << _sectionId
          << " time: "               << _time                 << "\n"
          << " energyDep: "          << _energyDep            << "\n"
          << " COG3Vector.u: "       << _cog3Vector.x()
          << " COG3Vector.v: "       << _cog3Vector.y()
          << " COG3Vector.w: "       << _cog3Vector.z()       <<"\n"
          << " size: "               << _caloCrystalHitsPtrVector.size()<<"\n";
    }

} 
