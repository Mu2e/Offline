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

#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/CrystalMapper.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <memory>




namespace mu2e {

   
    class Vane : public CaloSection {

 
      public:

          Vane(int vaneId, 
	       double rMin, 
	       int nCrystalX, 
	       int nCrystalY, 
	       const CLHEP::Hep3Vector &size,
	       double cellSize, 
	       double crystalHalfLength,
	       const CLHEP::Hep3Vector &vaneOriginToCrystalOrigin);



          std::vector<int> findLocalNeighbors(int crystalId, int level) const; 
          int              idxFromPosition(double x, double z)          const;

          void             print(std::ostream &os = std::cout)          const;



          virtual double innerRadius()    const  {return _rMin; }
          virtual double outerRadius()    const  {return _rMin+_nCrystalY*2.0*_cellSize;}
          int    crystalY(int id)         const  {return id/_nCrystalX;}
          int    crystalX(int id)         const  {return id%_nCrystalX;}
          double xActiveHalf()            const  {return _xActiveHalf; }
          double yActiveHalf()            const  {return _yActiveHalf; }



       private:
 
          void     fillCrystals(const CLHEP::Hep3Vector &crystalOriginInVane);

          double   _rMin;
          int      _nCrystalX;
          int      _nCrystalY;
          double   _cellSize;
	  double   _xActiveHalf;
	  double   _yActiveHalf;

    	  std::shared_ptr<CrystalMapper>  _crystalMap;

    };
}

#endif 
