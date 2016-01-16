// $Id: VaneMapper.cc,v 1.4 2013/07/25 23:56:46 echenard Exp $
// $Author: echenard $
// $Date: 2013/07/25 23:56:46 $
//
// Vane position map generator: 
//   tesselate a vane with crystals. For historical reasons, we start on the lower left corner of the vane
//


// C++ includes
#include <iostream>
#include <map>
#include <cmath>

// Mu2e includes
#include "CalorimeterGeom/inc/VaneMapper.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {


      VaneMapper::VaneMapper(int nCryX, int nCryY) : 
         _apexX(),
	 _apexY(),
	 _nCrystalX(nCryX),
	 _nCrystalY(nCryY)    
      {}
                       


      CLHEP::Hep2Vector VaneMapper::xyFromIndex(int thisIndex) const
      {        
         double x = (2*(thisIndex%_nCrystalX) - _nCrystalX+1);
         double y = (2*(thisIndex/_nCrystalX) - _nCrystalY+1);
	 return CLHEP::Hep2Vector(x,y);
      }
  
      int VaneMapper::indexFromXY(double x, double y) const
      {        
	  return  int(x/2.0 + 0.5*_nCrystalX) + int(y/2.0+0.5*_nCrystalY)*_nCrystalX ;
      }


      std::vector<int> VaneMapper::neighbors(int thisIndex, unsigned int level)  const
      {	 
          int X0     = thisIndex % _nCrystalX;
          int Y0     = thisIndex / _nCrystalX;
	  int ilevel = int(level);

          std::vector<int> list;
	  
	  //X0-level -> X0+level with Y=Y0+level
          for (int i=-ilevel;i<=ilevel;++i) {
            int X = X0+i, Y = Y0+ilevel;
            if (X < 0 || X > (_nCrystalX-1)  || Y > (_nCrystalY-1) ) continue;
            list.push_back(Y*_nCrystalX + X);
          }

          //Y0+ilevel-1 -> Y0-ilevel with X=X0+ilevel
          for (int i=ilevel-1;i>=-ilevel;--i) {
             int X = X0+ilevel, Y = Y0+i;
             if (Y < 0 || Y > (_nCrystalY-1)|| X > (_nCrystalX-1)  ) continue;
             list.push_back(Y*_nCrystalX + X); 
           }

           //X0+ilevel-1 -> X0-ilevel with Y=Y0-ilevel
           for (int i=ilevel-1;i>=-ilevel;--i) {
             int X = X0+i, Y = Y0-ilevel;
             if (X < 0 || X > (_nCrystalX-1) || Y<0) continue;
             list.push_back(Y*_nCrystalX + X); 
           }

           //Y0-ilevel+1 -> Y0-ilevel-1 with X=X0-ilevel
           for (int i=-ilevel+1;i<ilevel;++i) {
             int X = X0-ilevel, Y = Y0+i;
             if (Y < 0 || Y > (_nCrystalY-1) || X < 0) continue;
             list.push_back(Y*_nCrystalX + X); 
           }
           
	   return list;	   
      }



}


             
