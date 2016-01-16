//
// $Id: Vane.cc,v 1.4 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//
// Note: We define the xy position of the crystal relative to the disk local coordinate system
// x is the "long" side, y is the "short" side.


// C++ includes
#include <iostream>


// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/VaneMapper.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"



namespace mu2e {

      Vane::Vane(int vaneId, double rMin, int nCrystalX, int nCrystalY, CLHEP::Hep3Vector const& size, 
                 double cellSize, double crystalHalfLength, 
		 CLHEP::Hep3Vector const & vaneOriginToCrystalOrigin) : 
		 
        CaloSection(vaneId,size,vaneOriginToCrystalOrigin),
	_rMin(rMin),
	_nCrystalX(nCrystalX),
	_nCrystalY(nCrystalY), 
	_cellSize(cellSize),
	_xActiveHalf(nCrystalX*cellSize+1e-6),
	_yActiveHalf(nCrystalY*cellSize+1e-6)
	
       {
   	   _crystalMap = std::shared_ptr<CrystalMapper>(new VaneMapper(nCrystalX,nCrystalY));
           
	   fillCrystals(vaneOriginToCrystalOrigin); //see note in VaneCalorimeterMaker
       }

     
      void Vane::fillCrystals(CLHEP::Hep3Vector const &crystalOriginInVane)
      {         

          int nCrystals = _nCrystalX*_nCrystalY;	  
          for (int i=0; i< nCrystals; ++i)
	  {              
	      CLHEP::Hep2Vector xy  = _cellSize*_crystalMap->xyFromIndex(i);
              CLHEP::Hep3Vector posFF(xy.x(),xy.y(),0);
	      CLHEP::Hep3Vector pos = posFF + crystalOriginInVane;

              _crystalList.push_back( Crystal(i,_id, pos,posFF) );
          }

      }

      int Vane::idxFromPosition(double x, double y) const 
      {        
	 return  _crystalMap->indexFromXY(x/_cellSize,y/_cellSize);
      }



      //find the local index of the neighbors for a given level (level = number of "layers" away)
      std::vector<int> Vane::findLocalNeighbors(int crystalId, int level) const
      {
	 return _crystalMap->neighbors(crystalId,level);
      }


      void Vane::print(std::ostream &os) const 
      {        
	 os<<"Number of crystals X / Y= "<<_nCrystalX<<" / "<<_nCrystalY<<" / "<<_nCrystalX*_nCrystalY<<std::endl;
	 os<<"Radius In / Out "<<_rMin<<" / "<<_rMin+_nCrystalY*_cellSize<<std::endl;
      }


}

