#ifndef CalorimeterGeom_Vane_hh
#define CalorimeterGeom_Vane_hh
//
// $Id: Vane.hh,v 1.11 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
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

         double innerRadius(void) const              {return _rMin; }
         double outerRadius(void) const              {return _rMin+_nCrystalR*2.0*_cellSize;}

         int crystalR(int id) const               {return id/_nCrystalZ;}
         int crystalZ(int id) const               {return id%_nCrystalZ;}

         std::vector<int> neighbors(int crystalId, int level=1) const; 



       private:
 
         double               _rMin;
         int                  _nCrystalR;
         int                  _nCrystalZ;
         double               _cellSize;

         void fillCrystals(void);
         std::vector<int> findNeighbors(int crystalId, int level) const; 

    };
}

#endif /* CalorimeterGeom_Vane_hh */
