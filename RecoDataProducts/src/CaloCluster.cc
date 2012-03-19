//
// $Id: CaloCluster.cc,v 1.1 2012/03/19 19:35:42 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/03/19 19:35:42 $
//
// Original author G. Pezzullo
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib/exception.h"

// Mu2e includes
#include "RecoDataProducts/inc/CaloCluster.hh"


using namespace std;

namespace mu2e {

void CaloCluster::SetVaneId (int vane) {
        _vaneId = vane;
}
void CaloCluster::SetTime (float time) {
        _time = time;
}
void CaloCluster::SetCogRow ( int row) {
        _cogRow = row;
}
void CaloCluster::SetCogColumn ( int column) {
        _cogColumn = column;
}
void CaloCluster::SetEnergyDep ( float energyDep) {
        _energyDep = energyDep;
}
void CaloCluster::SetCog3Vector ( CLHEP::Hep3Vector cog3Vector) {
        _cog3Vector = cog3Vector;
}
void CaloCluster::SetCog3VectorError ( CLHEP::Hep3Vector cog3VectorErr) {
        _cog3VectorError = cog3VectorErr;
}
void CaloCluster::SetShowerDir( float dir){
        _showerDir = dir;
}

void CaloCluster::SetErrShowerDir(float errDir){
        _errShowerDir = errDir;
}

void CaloCluster::AddHit (CaloCrystalHitPtr &a) {

        _caloCrystalHitsPtrVector.push_back(a);
        //time*=(float)clusterSize;
        _time*=(float)_caloCrystalHitsPtrVector.size();
        //++clusterSize;
        _time += a->time();
        //time /= (float)clusterSize;
        _time /= (float)_caloCrystalHitsPtrVector.size();
        _energyDep += a->energyDep();
}

// Print the information found in this hit.
void CaloCluster::print( ostream& ost, bool doEndl ) const {

        ost << "CaloCluster :   "
                        << " vane: "          << _vaneId
                        << " time: "          << _time
                        << " COGrow: "        << _cogRow
                        << " energyDep: "     << _energyDep
                        << " COGcolumn: "     << _cogColumn
                        << " COG3Vector.u: "    << _cog3Vector.x()
                        << " COG3Vector.v: "    << _cog3Vector.y()
                        << " COG3Vector.w: "    << _cog3Vector.z()
                        << " COG3VectorError.u: "    << _cog3VectorError.x()
                        << " COG3VectorError.v: "    << _cog3VectorError.y()
                        << " COG3VectorError.w: "    << _cog3VectorError.z()
                        << " size: "          << _caloCrystalHitsPtrVector.size();

        if ( doEndl ){
                ost << endl;
        }

}

} // namespace mu2e
