#ifndef CalorimeterGeom_CaloSection_hh
#define CalorimeterGeom_CaloSection_hh

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
            _id(id), _origin(), _originLocal(), _rotation(), _size(),
            _crystalShift(crystalShift), _crystalList() 
         {} 

	 int id(void) const                                       {return _id;}
         int nCrystals(void) const                                {return _crystalList.size();}
         Crystal const& getCrystal(int i) const                   {return _crystalList.at(i);}

	 CLHEP::Hep3Vector  const& getOrigin(void) const          {return _origin;}
         CLHEP::Hep3Vector  const& getOriginLocal(void) const     {return _originLocal; }
	 CLHEP::HepRotation const& getRotation(void) const        {return _rotation;}
	 CLHEP::HepRotation const& getInverseRotation(void) const {return _rotationInverse;}
         CLHEP::Hep3Vector  const& size(void) const               {return _size; }

         void setOrigin(const CLHEP::Hep3Vector& orig)            {_origin = orig;}
         void setOriginLocal(const CLHEP::Hep3Vector& orig)       {_originLocal = orig;}
         void setRotation(const CLHEP::HepRotation& rot)          {_rotation = rot; _rotationInverse=rot.inverse();}
         void setSize(const CLHEP::Hep3Vector& size)              {_size = size;}

      protected:

	 int                  _id;
	 CLHEP::Hep3Vector    _origin;
         CLHEP::Hep3Vector    _originLocal;
	 CLHEP::HepRotation   _rotation;
	 CLHEP::HepRotation   _rotationInverse;
	 CLHEP::Hep3Vector    _size;

	 CLHEP::Hep3Vector    _crystalShift;
         std::vector<Crystal> _crystalList;



   };

}

#endif /* CalorimeterGeom_CaloSection_hh */
