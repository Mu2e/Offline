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
	     _frontFaceCenter(CLHEP::Hep3Vector(0,0,0)),
	     _backFaceCenter(CLHEP::Hep3Vector(0,0,0)), 
	     _innerRadius(0),
	     _outerRadius(0)
	   {}


	   int                        id()                     const {return _id;}
           
	   int                        nCrystals()              const {return _crystalList.size();}
           const Crystal&             crystal(int i)           const {return _crystalList.at(i);}
           Crystal&                   crystal(int i)                 {return _crystalList.at(i);}   


           const CLHEP::Hep3Vector&  size()                    const {return _size; }
	   const CLHEP::Hep3Vector&  origin()                  const {return _origin;}
           const CLHEP::Hep3Vector&  originLocal()             const {return _originLocal; }
	   const CLHEP::Hep3Vector&  originToCrystalOrigin()   const {return _originToCrystalOrigin;}
	   const CLHEP::HepRotation& rotation()                const {return _rotation;}
	   const CLHEP::HepRotation& inverseRotation()         const {return _inverseRotation;}
           const CLHEP::Hep3Vector&  frontFaceCenter()         const {return _frontFaceCenter; }
           const CLHEP::Hep3Vector&  backFaceCenter()          const {return _backFaceCenter; }
	   double innerEnvelopeR()                             const {return _innerRadius;}
	   double outerEnvelopeR()                             const {return _outerRadius;}
           
	   void setOrigin(const CLHEP::Hep3Vector &orig)          {_origin          = orig;}
           void setOriginLocal(const CLHEP::Hep3Vector &orig)     {_originLocal     = orig;}
           void setRotation(const CLHEP::HepRotation &rot)        {_rotation        = rot; _inverseRotation = rot.inverse();}
           void setFrontFaceCenter(const CLHEP::Hep3Vector &pos)  {_frontFaceCenter = pos;}
           void setBackFaceCenter(const CLHEP::Hep3Vector &pos)   {_backFaceCenter  = pos;}
           void setEnvelopeRad(double rin, double rout)           {_innerRadius=rin; _outerRadius=rout;}
           
           void print(std::ostream &os = std::cout) const;



       protected:

           std::vector<Crystal> _crystalList;

	   int                  _id;
	   CLHEP::Hep3Vector    _size;
	   CLHEP::Hep3Vector    _origin;
           CLHEP::Hep3Vector    _originLocal;
	   CLHEP::HepRotation   _rotation;
	   CLHEP::HepRotation   _inverseRotation;

	   CLHEP::Hep3Vector    _originToCrystalOrigin;
	   CLHEP::Hep3Vector    _frontFaceCenter;
	   CLHEP::Hep3Vector    _backFaceCenter;
	   
	   double               _innerRadius;
	   double               _outerRadius;
   };

}

#endif
