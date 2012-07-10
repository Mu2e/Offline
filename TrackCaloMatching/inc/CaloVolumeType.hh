//
//
//
// $Id: CaloVolumeType.hh,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef CALOVOLUMETYPE_HH
#define CALOVOLUMETYPE_HH


#include "BaBar/include/DetectorModel/DetVolumeType.hh"
#include "BaBar/include/DetectorModel/DetSurface.hh"
#include "MatEnv/MatDBInfo.hh"


namespace mu2e{
class CaloVolumeType : public DetVolumeType {
public :
        CaloVolumeType( const char* name, int id ):
                DetVolumeType( name, id ){
                _matDB =  new MatDBInfo();
                _material = _matDB->findDetMaterial("LYSO");
        }
        ~CaloVolumeType (){}



        //Setting parameters
        void AddSide(DetSurface *side);
        void AddSideCorner(SurfacePoint *sideCorner, unsigned int sideId);

        //  Determine whether a given (local) coordinate lies 'within' the object.
        bool physicalMaterial(const TypeCoord*) const {return true;} // physical size
        //  Describe the material at a given coordinate
        const DetMaterial& material(const TypeCoord*) const {return *_material;}


protected:

        // Helper functions
        bool insideLine( const SurfacePoint& thisPoint,
                        const SurfacePoint& p1,
                        const SurfacePoint& p2 ) const ;
private :

        const DetMaterial *_material;
        MatDBInfo* _matDB;

};

}

#endif
