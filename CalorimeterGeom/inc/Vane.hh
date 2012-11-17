#ifndef CalorimeterGeom_Vane_hh
#define CalorimeterGeom_Vane_hh

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
         double outerRadius(void) const              {return _rMin+_nCrystalR*_cellSize;}

         int getCrystalR(int id) const               {return id/_nCrystalZ;}
         int getCrystalZ(int id) const               {return id%_nCrystalZ;}

         std::vector<int> getNeighbors(int crystalId, int level=1) const; 



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
