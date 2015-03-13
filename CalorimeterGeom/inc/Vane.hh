#ifndef CalorimeterGeom_Vane_hh
#define CalorimeterGeom_Vane_hh
//
// $Id: Vane.hh,v 1.13 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Hold information about a vane in the calorimter.
//
// Original author R. Kutschke, Modified B. Echenard
//

//C++ includes
#include <vector>

//mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Crystal.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e {

   
    class Vane : public CaloSection {

 
      public:

          Vane(int vaneId, double rMin, int nCrystalX, int nCrystalY, CLHEP::Hep3Vector const& size,
	       double cellSize, double crystalHalfLength,
	       CLHEP::Hep3Vector const& vaneOriginToCrystalOrigin);



          std::vector<int> findLocalNeighbors(int crystalId, int level) const; 
          int              idxFromPosition(double x, double z)          const;

          void             print()                                      const;



          double innerRadius()    const    {return _rMin; }
          double outerRadius()    const    {return _rMin+_nCrystalY*2.0*_cellSize;}
          int    crystalY(int id) const    {return id/_nCrystalX;}
          int    crystalX(int id) const    {return id%_nCrystalX;}
          double xActiveHalf()    const    {return _xActiveHalf; }
          double yActiveHalf()    const    {return _yActiveHalf; }



       private:
 
          void     fillCrystals(CLHEP::Hep3Vector const& crystalOriginInVane);

          double   _rMin;
          int      _nCrystalX;
          int      _nCrystalY;
          double   _cellSize;
	  
	  double   _xActiveHalf;
	  double   _yActiveHalf;


    };
}

#endif /* CalorimeterGeom_Vane_hh */
