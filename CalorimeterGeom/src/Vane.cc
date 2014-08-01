//
// $Id: Vane.cc,v 1.4 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

// Note: crystalShift indicates the shift of the crystal center from the center of 
// the cell in the calo due to the readout (or other things).


// C++ includes
#include <iostream>


// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Vane.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e {

      Vane::Vane(int id, double rMin, int nCrystalR, int nCrystalZ, double cellSize, CLHEP::Hep3Vector crystalShift) : 
        CaloSection(id, crystalShift),
        _rMin(rMin), _nCrystalR(nCrystalR), _nCrystalZ(nCrystalZ),_cellSize(cellSize)
      {
        fillCrystals();
      }

     
      void Vane::fillCrystals(void)
      {         

          int nCrystals = _nCrystalZ*_nCrystalR;	  
          for (int i=0; i< nCrystals; ++i)
	  {              
              double y = (2*(i/_nCrystalZ) - _nCrystalR+1)*_cellSize;
              double z = (2*(i%_nCrystalZ) - _nCrystalZ+1)*_cellSize;
              CLHEP::Hep3Vector pos(0,y,z);
              pos += _crystalShift;  //see note 
              _crystalList.push_back( Crystal(i,pos) );
          }

      }

      int Vane::idxFromPosition(double y, double z) const 
      {        
	return  int(y/2.0/_cellSize+0.5*_nCrystalR)*_nCrystalZ + int(z/2.0/_cellSize + 0.5*_nCrystalZ);
      }


      //find the local index of the neighbors for a given level (level = number of "layers" away)
      std::vector<int> Vane::findLocalNeighbors(int crystalId, int level) const
      {
          int Z0 = crystalId % _nCrystalZ;
          int R0 = crystalId / _nCrystalZ;

          std::vector<int> list;
	  
	  //Z0-level -> Z0+level with R=R0+level
          for (int i=-level;i<=level;++i) {
            int Z = Z0+i, R = R0+level;
            if (Z < 0 || Z > (_nCrystalZ-1)  || R > (_nCrystalR-1) ) continue;
            list.push_back(R*_nCrystalZ + Z);
          }

          //R0+level-1 -> R0-level with Z=Z0+level
          for (int i=level-1;i>=-level;--i) {
             int Z = Z0+level, R = R0+i;
             if (R < 0 || R > (_nCrystalR-1)|| Z > (_nCrystalZ-1)  ) continue;
             list.push_back(R*_nCrystalZ + Z); 
           }

           //Z0+level-1 -> Z0-level with R=R0-level
           for (int i=level-1;i>=-level;--i) {
             int Z = Z0+i, R = R0-level;
             if (Z < 0 || Z > (_nCrystalZ-1) || R<0) continue;
             list.push_back(R*_nCrystalZ + Z); 
           }

           //R0-level+1 -> R0-level-1 with Z=Z0-level
           for (int i=-level+1;i<level;++i) {
             int Z = Z0-level, R = R0+i;
             if (R < 0 || R > (_nCrystalR-1) || Z < 0) continue;
             list.push_back(R*_nCrystalZ + Z); 
           }
           
	   return list;	   
      }


}

