#ifndef CalorimeterGeom_Disk_hh
#define CalorimeterGeom_Disk_hh
// $Id: Disk.hh,v 1.4 2013/03/08 01:22:31 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:31 $
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

          std::vector<int> neighbors(int crystalId, int level=1) const; 
          int              idxFromPosition(double x, double y) const;
          double           estimateEmptySpace() const;

          double innerRadius() const   {return _radiusIn;}
          double outerRadius() const   {return _radiusOut;}
          double thickness()   const   {return _thickness;}


       private:

          void             fillCrystals();
          bool             isInsideDisk(double x, double y) const;
	  double           calcDistToSide(CLHEP::Hep2Vector& a, CLHEP::Hep2Vector& b) const;
          std::vector<int> findNeighbors(int crystalId, int level=1) const; 


	  double           _radiusIn;
	  double           _radiusOut;
          double           _thickness;
	  double           _cellSize;	 
          HexMap           _hexMap;
          std::vector<int> _mapToCrystal;
	  std::vector<int> _crystalToMap;

    };
}

#endif /* CalorimeterGeom_Disk_hh */
