#ifndef CalorimeterGeom_Disk_hh
#define CalorimeterGeom_Disk_hh
// $Id: Disk.hh,v 1.3 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Hold information about a disk in the calorimter.
//
// Original author B Echenard 
//

// C++ includes
#include <vector>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/HexMap.hh"
#include "CalorimeterGeom/inc/Crystal.hh"

//CLHEP includess
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {


    class Disk : public CaloSection {


       public:

	  Disk(int id, double rin, double rout, double thickness, double cellSize, CLHEP::Hep3Vector crystalOffset); 

          double innerRadius(void) const           {return _radiusIn;}
          double outerRadius(void) const           {return _radiusOut;}
          double thickness(void) const             {return _thickness;}

          std::vector<int> neighbors(int crystalId, int level=1) const; 
          double EstimateEmptySpace(void);



       private:

	  double               _radiusIn;
	  double               _radiusOut;
          double               _thickness;
	  double               _cellSize;
	 
          HexMap               _HexMap;
          std::vector<int>     _mapToCrystal;
	  std::vector<int>     _CrystalToMap;

          void   fillCrystals(void);
          bool   isInsideDisk(double x, double y) const;
	  double calcDistToSide(CLHEP::Hep2Vector& a, CLHEP::Hep2Vector& b) const;
          std::vector<int> findNeighbors(int crystalId, int level=1) const; 
    };
}

#endif /* CalorimeterGeom_Disk_hh */
