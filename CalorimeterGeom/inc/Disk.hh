#ifndef CalorimeterGeom_Disk_hh
#define CalorimeterGeom_Disk_hh

//
// Hold information about a disk in the calorimter.
//
// Original author B Echenard 
//

// C++ includes
#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

// Mu2e includes
#include "CalorimeterGeom/inc/HexPositionMap.hh"



namespace mu2e {


   class Disk{


      public:


	 Disk(double rin, double rout, double thick, double cellSize, int id=-1) : 
	 _radiusIn(rin),_radiusOut(rout),_thickness(thick),_cellSize(cellSize),_id(id),_origin(),
	 _rotation(),_crystalMap(cellSize, rin, rout)
	 {            
	 }

	 int    id(void) const                                 {return _id;}
	 double getRin(void) const                             {return _radiusIn;}
	 double getRout(void) const                            {return _radiusOut;}
    
	 CLHEP::Hep3Vector const& getOrigin(void) const        {return _origin;}
	 CLHEP::HepRotation const& getRotation(void) const     {return _rotation;}
	 void setOrigin(const CLHEP::Hep3Vector& orig)         {_origin   = orig;}
	 void setRotation(const CLHEP::HepRotation& rot)       {_rotation = rot;}
         
	 HexPositionMap const& getCrystalMap(void) const       {return _crystalMap;} 
         

      private:

	 double _radiusIn;
	 double _radiusOut;
	 double _thickness;
	 double _cellSize;
	 int    _id;

	 CLHEP::Hep3Vector  _origin;
	 CLHEP::HepRotation _rotation;
  
	 HexPositionMap _crystalMap;
  };


} //namespace mu2e


#endif /* CalorimeterGeom_Disk_hh */
