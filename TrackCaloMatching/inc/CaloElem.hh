//
//
//
// $Id: CaloElem.hh,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef CALOELEM_HH
#define CALOELEM_HH


//#include "BaBar/include/DetectorModel/DetVolumeType.hh"
#include "BaBar/include/CLHEP/Geometry/HepPoint.h"
//#include "BaBar/include/CLHEP/Geometry/Transformation.h"
#include "BaBar/BaBar/include/Constants.hh"
#include "BaBar/include/DetectorModel/DetElem.hh"
#include "TrackCaloMatching/inc/Calorimeter4VanesGeom.hh"
#include "TrkBase/TrkRep.hh"
#include <cmath>
//#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
class HepPoint;
//class HepTransformation;

namespace mu2e{
class CaloElem : public DetElem {
public :

        //Constructors
        CaloElem(Calorimeter4VanesGeom& caloVanes, TrkRep const* trep);
        CaloElem(Calorimeter4VanesGeom& caloVanes, TrkRep const* trep, const DetType*,const char*,int);

        int intersect(const Trajectory*,DetIntersection&) const;



        HepPoint coordToPoint( const TypeCoord* aCoord ) const;

private :
        int _nVanes;
        Calorimeter4VanesGeom _CaloVanes;
        TrkRep const* _trep;
};
}


#endif
