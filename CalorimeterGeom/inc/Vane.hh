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

          Vane(int id, double rMin, int nCrystalR, int nCrystalZ, double cellSize, CLHEP::Hep3Vector crystalOffset);

          std::vector<int> findLocalNeighbors(int crystalId, int level) const; 
          int              idxFromPosition(double y, double z)          const;


          double innerRadius() const    {return _rMin; }
          double outerRadius() const    {return _rMin+_nCrystalR*2.0*_cellSize;}
          int crystalR(int id) const    {return id/_nCrystalZ;}
          int crystalZ(int id) const    {return id%_nCrystalZ;}


       private:
 
          void fillCrystals();

          double   _rMin;
          int      _nCrystalR;
          int      _nCrystalZ;
          double   _cellSize;


    };
}

#endif /* CalorimeterGeom_Vane_hh */
