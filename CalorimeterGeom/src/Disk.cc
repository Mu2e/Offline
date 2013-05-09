//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

//Notes: HexMap tesselates a plane with hexagons, we need a ring
// _mapToCrystal and _crystalToMap get the crystal number conversion 
// between the map and the ring (crystal number = position in crystalList)
//
// CrystalShift indicates the shift of the crystal center from the center of 
// the cell in the calo due to the readout and other things. This is usually disk/vane dependent


// C++ includes
#include <iostream>
#include <cmath>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/HexMap.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

      Disk::Disk(int id, double rin, double rout, double cellSize, CLHEP::Hep3Vector crystalShift) : 
        CaloSection(id, crystalShift),_radiusIn(rin), _radiusOut(rout),
	_cellSize(cellSize), _hexMap(),_mapToCrystal(),_crystalToMap()
      { 
        fillCrystals(); 
      }

     

      // take the crystals from the HexMap, and keep only those who are in the annulus
      void Disk::fillCrystals(void)
      {   

	    int nRings    = int(1.5*_radiusOut/_cellSize)+1;
	    int nHexagons = 1 + 3*nRings*(nRings-1);

	    _mapToCrystal.reserve(nHexagons);
	    _crystalToMap.reserve(nHexagons);

            int nCrystal(0);
 	    for (int i=0;i<nHexagons;++i) {

                CLHEP::Hep2Vector xy = _hexMap.xyPosition(i);

 	        if ( !isInsideDisk(_cellSize*xy.x(),_cellSize*xy.y()) ) {_mapToCrystal.push_back(-1); continue;}

                _crystalToMap.push_back(i);
     	        _mapToCrystal.push_back(nCrystal);

                CLHEP::Hep3Vector pos(_cellSize*xy.x(),_cellSize*xy.y(),0);
                pos += _crystalShift;

	        _crystalList.push_back( Crystal(nCrystal,pos) );	        
                ++nCrystal;
	    }


	    // precalculate nearest neighbours list for each crystal
	    for (unsigned int i=0;i<_crystalList.size();++i){
	 	  _crystalList[i].setNearestNeighbours(findNeighbors(i,1));
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
        return _mapToCrystal[mapIdx];
      }





      std::vector<int> Disk::neighbors(int crystalId, int level) const
      {
	   if (level==1) return _crystalList.at(crystalId).nearestNeighbours();
	   return findNeighbors(crystalId,level);
      }
      
      std::vector<int> Disk::findNeighbors(int crystalId, int level) const
      {

	   std::vector<int> list; 
	   list.reserve(20);

	   std::vector<int> temp( _hexMap.neighbors(_crystalToMap[crystalId],level) );
	   for (unsigned int i=0;i<temp.size();++i) {
	     if (_mapToCrystal[temp[i]] >-1) list.push_back(_mapToCrystal[temp[i]]);
	   }

	   return list;
      }


      double Disk::estimateEmptySpace(void) const
      {
	    double sum(0),dx(0.02),dy(0.02);
	    for (double x=0;x<=1.5*_radiusIn;x+=dx){

              double y0 = (x<_radiusIn) ? sqrt(_radiusIn*_radiusIn-x*x) : 0;
	      for (double y=y0;y<=y0+5*_cellSize;y+=dy){
		 int mapIdx = _hexMap.indexFromXY(x/_cellSize,y/_cellSize);
		 int iCry   = _mapToCrystal[mapIdx];
		 if (iCry==-1) sum+=dx*dy;		 
	      }  
	    }
            return 4.0*sum;
      }


}

