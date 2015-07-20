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
// x is the "long" side, y is the short side.


// C++ includes
#include <iostream>


// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Vane.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e {

      Vane::Vane(int vaneId, double rMin, int nCrystalX, int nCrystalY, CLHEP::Hep3Vector const& size, 
                 double cellSize, double crystalHalfLength, 
		 CLHEP::Hep3Vector const& vaneOriginToCrystalOrigin) : 
		 
        CaloSection(vaneId,size,vaneOriginToCrystalOrigin),
	_rMin(rMin),
	_nCrystalX(nCrystalX),
	_nCrystalY(nCrystalY), 
	_cellSize(cellSize)
       {
           fillCrystals(vaneOriginToCrystalOrigin); //see note in VaneCalorimeterMaker

 	   _xActiveHalf = nCrystalX*cellSize+1e-6;   
	   _yActiveHalf = nCrystalY*cellSize+1e-6;   
      }

     
      void Vane::fillCrystals(CLHEP::Hep3Vector const& crystalOriginInVane)
      {         

          int nCrystals = _nCrystalX*_nCrystalY;	  
          for (int i=0; i< nCrystals; ++i)
	  {              
              double x = (2*(i%_nCrystalX) - _nCrystalX+1)*_cellSize;
              double y = (2*(i/_nCrystalX) - _nCrystalY+1)*_cellSize;

              CLHEP::Hep3Vector pos(x,y,0);
              pos += crystalOriginInVane; 
              _crystalList.push_back( Crystal(i,_id, pos) );
          }

      }

      int Vane::idxFromPosition(double x, double y) const 
      {        
	return  int(x/2.0/_cellSize + 0.5*_nCrystalX) + int(y/2.0/_cellSize+0.5*_nCrystalY)*_nCrystalX ;
      }


      //find the local index of the neighbors for a given level (level = number of "layers" away)
      std::vector<int> Vane::findLocalNeighbors(int crystalId, int level) const
      {
          int X0 = crystalId % _nCrystalX;
          int Y0 = crystalId / _nCrystalX;

          std::vector<int> list;
	  
	  //X0-level -> X0+level with Y=Y0+level
          for (int i=-level;i<=level;++i) {
            int X = X0+i, Y = Y0+level;
            if (X < 0 || X > (_nCrystalX-1)  || Y > (_nCrystalY-1) ) continue;
            list.push_back(Y*_nCrystalX + X);
          }

          //Y0+level-1 -> Y0-level with X=X0+level
          for (int i=level-1;i>=-level;--i) {
             int X = X0+level, Y = Y0+i;
             if (Y < 0 || Y > (_nCrystalY-1)|| X > (_nCrystalX-1)  ) continue;
             list.push_back(Y*_nCrystalX + X); 
           }

           //X0+level-1 -> X0-level with Y=Y0-level
           for (int i=level-1;i>=-level;--i) {
             int X = X0+i, Y = Y0-level;
             if (X < 0 || X > (_nCrystalX-1) || Y<0) continue;
             list.push_back(Y*_nCrystalX + X); 
           }

           //Y0-level+1 -> Y0-level-1 with X=X0-level
           for (int i=-level+1;i<level;++i) {
             int X = X0-level, Y = Y0+i;
             if (Y < 0 || Y > (_nCrystalY-1) || X < 0) continue;
             list.push_back(Y*_nCrystalX + X); 
           }
           
	   return list;	   
      }

      void Vane::print() const 
      {        
	std::cout<<"Number of crystals X / Y= "<<_nCrystalX<<" / "<<_nCrystalY<<" / "<<_nCrystalX*_nCrystalY<<std::endl;
	std::cout<<"Radius In / Out "<<_rMin<<" / "<<_rMin+_nCrystalY*_cellSize<<std::endl;
      }


}

