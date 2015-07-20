#ifndef CalorimeterGeom_CaloSection_hh
#define CalorimeterGeom_CaloSection_hh
// $Id: CaloSection.hh,v 1.5 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
//
// Hold information about a CaloCaloSection in the calorimter.
//
// Original author B Echenard 
//

// C++ includes
#include <vector>

// Mu2e includes
#include "CalorimeterGeom/inc/Crystal.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"


namespace mu2e {


   class CaloSection {

       public:

	   CaloSection(int id, CLHEP::Hep3Vector const& size, CLHEP::Hep3Vector const& originToCrystalOrigin) : 
	     _crystalList(), 
	     _id(id), 
	     _size(size),
	     _origin(CLHEP::Hep3Vector(0,0,0)), 
	     _originLocal(CLHEP::Hep3Vector(0,0,0)),
	     _rotation(CLHEP::HepRotation::IDENTITY),
	     _inverseRotation(CLHEP::HepRotation::IDENTITY),
	     _originToCrystalOrigin(originToCrystalOrigin), 
	     _zDownInTracker(0),
	     _zUpInTracker(0),
	     _rInTracker(0),
	     _rOutTracker(0)
	   {}


	   int                        id()                       const {return _id;}
           
	   int                        nCrystals()                const {return _crystalList.size();}
           Crystal             const& crystal(int i)             const {return _crystalList.at(i);}
           Crystal&                   crystal(int i)                   {return _crystalList.at(i);}


           CLHEP::Hep3Vector  const& size()                      const {return _size; }
	   CLHEP::Hep3Vector  const& origin()                    const {return _origin;}
           CLHEP::Hep3Vector  const& originLocal()               const {return _originLocal; }
	   CLHEP::Hep3Vector  const& originToCrystalOrigin()     const {return _originToCrystalOrigin;}
	   CLHEP::HepRotation const& rotation()                  const {return _rotation;}
	   CLHEP::HepRotation const& inverseRotation()           const {return _inverseRotation;}

	   double zUpInTracker()                                 const {return _zUpInTracker;}
	   double zDownInTracker()                               const {return _zDownInTracker;}
	   double rInTracker()                                   const {return _rInTracker;}
	   double rOutTracker()                                  const {return _rOutTracker;}

           
	   void setOrigin(          const CLHEP::Hep3Vector& orig)  {_origin = orig;      }
           void setOriginLocal(     const CLHEP::Hep3Vector& orig)  {_originLocal = orig; }
           void setRotation(        const CLHEP::HepRotation& rot)  {_rotation = rot; _inverseRotation = rot.inverse();}

           void setBoundsInTracker(CLHEP::Hep3Vector  const& trackerOffset, double z0, double z1, double r0, double r1);
           void print() const;



       protected:

           std::vector<Crystal> _crystalList;

	   int                  _id;
	   CLHEP::Hep3Vector    _size;
	   CLHEP::Hep3Vector    _origin;
           CLHEP::Hep3Vector    _originLocal;
	   CLHEP::HepRotation   _rotation;
	   CLHEP::HepRotation   _inverseRotation;

	   CLHEP::Hep3Vector    _originToCrystalOrigin;
	   
	   double               _zDownInTracker;
	   double               _zUpInTracker;
	   double               _rInTracker;
	   double               _rOutTracker;
   };

}

#endif
