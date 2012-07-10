//
//
//
// $Id: CaloVolumeElem.hh,v 1.1 2012/07/10 00:02:19 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 00:02:19 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef CALOVOLUMEELEM_HH_
#define CALOVOLUMEELEM_HH_

#include "BaBar/include/DetectorModel/DetVolumeElem.hh"
#include "BaBar/include/DetectorModel/DetVolumeType.hh"
#include "TrackCaloMatching/inc/CaloVolumeType.hh"

class DetSurface;
class DetVolumeType;

class DetVolumeElem;
class Trajectory;
class CaloVolumeType;


namespace mu2e{
class CaloVolumeElem : public DetVolumeElem{
public :
        CaloVolumeElem( DetVolumeType* itsType, const char* name, int id,
                        const HepTransformation& theAlignment ):
                DetVolumeElem( itsType, name, id, theAlignment ){
                _sides = const_cast< std::vector< DetSurface* >* >(sides());
        }
        ~CaloVolumeElem (){}

        //Modifiers
        void createCache();
        void AddSurface(DetSurface*a)const
        {
               //sides()->push_back(a);
                _sides->push_back(a);
        }

private :
        std::vector< DetSurface* >* _sides;
};

}

#endif
