//
//
//
// $Id: CaloVolumeType.cc,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo & G. Tassielli
//

#include "TrackCaloMatching/inc/CaloVolumeType.hh"
#include <iostream>

namespace mu2e{
void CaloVolumeType::AddSide(DetSurface *side){
        mySides()->push_back(side);
}

void CaloVolumeType::AddSideCorner(SurfacePoint *sideCorner, unsigned int sideId){
        if(sideId < mySideCorners()->size()){
                mySideCorners()->at(sideId).push_back(sideCorner);
        }else {
                std::vector< SurfacePoint* > vec;
                mySideCorners()->push_back(vec);
                mySideCorners()->at(sideId).push_back(sideCorner);
        }
}

bool CaloVolumeType::insideLine( const SurfacePoint& thisPoint, const SurfacePoint& p1, const SurfacePoint& p2 ) const {

        //that comparison works only if the center of the surface is also the origin of the local surface frame
        if( (std::fabs( thisPoint[0]) <= std::fabs(p1[0]) ) && (std::fabs( thisPoint[1]) <= std::fabs(p1[1]) ) ){
                return true;
        }else {
              return false;
        }
}

//const DetMaterial& material(const TypeCoord*) const


}


