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


   class CaloSection{

       public:

	   CaloSection(int id, CLHEP::Hep3Vector crystalShift) : 
              _crystalShift(crystalShift), _crystalList(), _id(id), _origin(), _originLocal(), _rotation(), _size()
           {}

	   int id() const                                           {return _id;}
           
	   int nCrystals() const                                    {return _crystalList.size();}
           Crystal const& crystal(int i) const                      {return _crystalList.at(i);}
           Crystal&       crystal(int i)                            {return _crystalList.at(i);}

	   CLHEP::Hep3Vector  const& origin()          const        {return _origin;}
           CLHEP::Hep3Vector  const& originLocal()     const        {return _originLocal; }
	   CLHEP::HepRotation const& rotation()        const        {return _rotation;}
	   CLHEP::HepRotation const& inverseRotation() const        {return _rotationInverse;}
           CLHEP::Hep3Vector  const& size()            const        {return _size; }
           CLHEP::Hep3Vector  const& crystalShift()    const        {return _crystalShift; }

           void setOrigin(const CLHEP::Hep3Vector& orig)            {_origin = orig;}
           void setOriginLocal(const CLHEP::Hep3Vector& orig)       {_originLocal = orig;}
           void setRotation(const CLHEP::HepRotation& rot)          {_rotation = rot; _rotationInverse = rot.inverse();}
           void setSize(const CLHEP::Hep3Vector& size)              {_size = size;}
           void setCrystalShift(const CLHEP::Hep3Vector& shift)     {_crystalShift = shift;}



       protected:

	   CLHEP::Hep3Vector    _crystalShift;
           std::vector<Crystal> _crystalList;

	   int                  _id;
	   CLHEP::Hep3Vector    _origin;
           CLHEP::Hep3Vector    _originLocal;
	   CLHEP::HepRotation   _rotation;
	   CLHEP::HepRotation   _rotationInverse;
	   CLHEP::Hep3Vector    _size;




   };

}

#endif /* CalorimeterGeom_CaloSection_hh */
