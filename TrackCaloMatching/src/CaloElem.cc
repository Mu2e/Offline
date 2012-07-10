//
// $Id: CaloElem.cc,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo & G. Tassielli
//

#include "TrackCaloMatching/inc/CaloElem.hh"
//#include "art/Utilities/Exception.h"
//#include <cassert>

using namespace std;

namespace mu2e{

CaloElem::CaloElem(Calorimeter4VanesGeom& caloVanes, TrkRep const* trep):DetElem(){
        _CaloVanes = caloVanes;
        _nVanes    = caloVanes.nVanes();
        _trep      = trep;
}

CaloElem::CaloElem(Calorimeter4VanesGeom& caloVanes, TrkRep const* trep,
                const DetType* type,const char* name,int id):DetElem(type, name, id){
        _CaloVanes = caloVanes;
        _nVanes    = caloVanes.nVanes();
        _trep      = trep;
}

int CaloElem::intersect(const Trajectory* traj, DetIntersection& dinter) const{
       // HelixTraj trkHel(_trep->helix(0).params(),_trep->helix(0).covariance());

        //double angle = Constants::pi*0.5 + trkHel.phi0();

        //double circleRadius = 1.0/trkHel.omega();
        //double centerCircleX = trkHel.d0() + circleRadius;
        //double centerCircleY = centerCircleX*sin(angle);
        //centerCircleX *= cos(angle);

        //double lowrange = trkHel.zFlight(_CaloVanes.ZfrontFaceCalo() ), highrange = trkHel.zFlight(_CaloVanes.ZbackFaceCalo() );
        int res0 = -1;
        //Length length[_nVanes];
        //_CaloVanes.caloExtrapol(_trep, lowrange, highrange, trkHel, res0, dinter, length);

        return res0;
}

HepPoint CaloElem::coordToPoint( const TypeCoord* aCoord ) const{
        HepPoint aPoint( (*aCoord)[0], (*aCoord)[1], (*aCoord)[2] );
        return transform().transFrom( aPoint );
}

}
