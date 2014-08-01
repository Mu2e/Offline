//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

//Notes: HexMap tesselates a plane with hexagons, we need a ring
// _posUtil handle number conversion between the map and the disk numbering scheme
//
// CrystalShift indicates the shift of the crystal center from the center of 
// the cell in the calo due to the readout and other things. This is usually disk/vane dependent


// C++ includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/HexMap.hh"
#include "CalorimeterGeom/inc/DiskCrystalPosUtil.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

      Disk::Disk(int id, double rin, double rout, double cellSize, CLHEP::Hep3Vector crystalShift) : 
        CaloSection(id, crystalShift),_radiusIn(rin), _radiusOut(rout),
	_cellSize(cellSize), _hexMap(),_posUtil()
      { 
        fillCrystals(); 
      }

     

      // take the crystals from the HexMap, and keep only those who are in the annulus
      void Disk::fillCrystals(void)
      {   

	    int nRings    = int(1.5*_radiusOut/_cellSize)+1;
	    int nHexagons = 1 + 3*nRings*(nRings-1);

            int nCrystal(0);
 	    for (int i=0;i<nHexagons;++i)
	    {
                CLHEP::Hep2Vector xy = _hexMap.xyFromIndex(i);

 	        if ( !isInsideDisk(_cellSize*xy.x(),_cellSize*xy.y()) ) {_posUtil.Fill(-1); continue;}

                CLHEP::Hep3Vector pos(_cellSize*xy.x(),_cellSize*xy.y(),0);
                pos += _crystalShift;

		_posUtil.Fill(i,nCrystal,_hexMap.l(i),_hexMap.k(i));
	        _crystalList.push_back( Crystal(nCrystal,pos) );
		
                ++nCrystal;
	    }
      }


      bool Disk::isInsideDisk(double x, double y) const
      {    	 
         // x,y coordinates of apexes of unit hexagon (position 7 is first apex again to close the loop)
         double apexX[7]={-0.2886751,+0.2886751,+0.5773502,+0.2886751,-0.2886751,-0.5773502,-0.2886751};
         double apexY[7]={-0.5,-0.5,0,0.5,0.5,0,-0.5};
         
         for (int i=0;i<6;++i) {  

             CLHEP::Hep2Vector p1(x+_cellSize*apexX[i],   y+_cellSize*apexY[i]);
             CLHEP::Hep2Vector p2(x+_cellSize*apexX[i+1], y+_cellSize*apexY[i+1]);

      	     //check distances. Note that the farthest distance is always at an apex in our case 
             if (calcDistToSide(p1,p2) < _radiusIn) return false;
             if (p1.mag() > _radiusOut)             return false;       
         }

         return true;
      }

      double Disk::calcDistToSide(CLHEP::Hep2Vector& a, CLHEP::Hep2Vector& b) const
      {
	  double t = -3.0*(b-a).dot(a);
	  if (t<0) return a.mag();
	  if (t>1) return b.mag();
          return (a+t*(b-a)).mag();
      }




      int Disk::idxFromPosition(double x, double y) const 
      {
        int mapIdx = _hexMap.indexFromXY(x/_cellSize,y/_cellSize);
        return _posUtil.mapToCrystal(mapIdx);
      }





      //find the local index of the neighbors for a given level (level = number of rings away)
      std::vector<int> Disk::findLocalNeighbors(int crystalId, int level) const
      {
	   std::vector<int> list; 
	   std::vector<int> temp(_hexMap.neighbors(_posUtil.crystalToMap(crystalId),level));

	   for (unsigned int i=0;i<temp.size();++i)
	     if (_posUtil.mapToCrystal(temp[i])>-1) list.push_back(_posUtil.mapToCrystal(temp[i]));
	   
	   return list;
      }

      
      
      //Slightly inefficient but robust integration. Divide the area between the disk and 
      //the first few crystals into tiny squares, and sum squares outside the crystals
      double Disk::estimateEmptySpace(void) const
      {
	    double sum(0),delta(0.1);
	    double RadiusMax=_radiusIn + 0.2*(_radiusOut-_radiusIn);
	    
	    for (double x=delta/2.0; x <= RadiusMax; x+=delta)
	    {
              double y0 = (x < _radiusIn) ? sqrt(_radiusIn*_radiusIn - x*x) + delta/2.0  : delta/2.0;
	      double ymax = sqrt(RadiusMax*RadiusMax-x*x);
	      for (double y=y0; y <= ymax; y+=delta)
	      {
		 int mapIdx = _hexMap.indexFromXY( x/_cellSize,y/_cellSize);
		 int iCry   = _posUtil.mapToCrystal(mapIdx);
		 if (iCry==-1) sum+=delta*delta;		 
	      }  
	    }
	    return 4.0*sum;
      }

    
}

