//
// $Id: CaloCluster.cc,v 1.3 2012/09/06 19:58:26 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/09/06 19:58:26 $
//
// Original author G. Pezzullo
//

// C++ includes
#include <ostream>

// Mu2e includes
#include "RecoDataProducts/inc/CaloCluster.hh"

using namespace std;

namespace mu2e {

  void CaloCluster::AddHit (CaloCrystalHitPtr &a) {

    _caloCrystalHitsPtrVector.push_back(a);

    _time *=_energyDep;
    _time += (a->time())*(a->energyDep());

    _energyDep += a->energyDep();
    _time /= _energyDep;
  }

  // Print the information found in this hit.
  void CaloCluster::print( ostream& ost, bool doEndl ) const {

    ost << "CaloCluster :   "
        << " vane: "               << _vaneId
        << " time: "               << _time                 << "\n"
        << " energyDep: "          << _energyDep            << "\n"
        << " COGrow: "             << _cogRow
        << " COGcolumn: "          << _cogColumn            <<"\n"
        << " COG3Vector.u: "       << _cog3Vector.x()
        << " COG3Vector.v: "       << _cog3Vector.y()
        << " COG3Vector.w: "       << _cog3Vector.z()       <<"\n"
        << " COG3VectorError.u: "  << _cog3VectorError.x()
        << " COG3VectorError.v: "  << _cog3VectorError.y()
        << " COG3VectorError.w: "  << _cog3VectorError.z()
        << " size: "               << _caloCrystalHitsPtrVector.size();

    if ( doEndl ){
      ost << endl;
    }

  }

} // namespace mu2e
