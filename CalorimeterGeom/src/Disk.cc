//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

// Notes: CrystalMap tesselates a plane with hexagons or squares, then selects crystals inside the ring



#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/HexMapper.hh"
#include "CalorimeterGeom/inc/SquareMapper.hh"
#include "CalorimeterGeom/inc/SquareShiftMapper.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"




namespace mu2e {


      Disk::Disk(int diskId, double rin, double rout, CLHEP::Hep3Vector const& size, 
                 double cellSize, int crystalNedges, bool shiftCrystal, double crystalHalfLength, 
		 CLHEP::Hep3Vector const& diskOriginToCrystalOrigin) :         
	 CaloSection(diskId,size,diskOriginToCrystalOrigin),
	 _radiusIn(rin),
	 _radiusOut(rout),
	 _cellSize(cellSize),
	 _mapToCrystal(),
	 _crystalToMap(),
         _crateDeltaZ(0)
      { 
	   if (crystalNedges==6) _crystalMap = std::shared_ptr<CrystalMapper>(new HexMapper());  
	   else
	   {
	      if (shiftCrystal) _crystalMap = std::shared_ptr<CrystalMapper>(new SquareShiftMapper());
	      else              _crystalMap = std::shared_ptr<CrystalMapper>(new SquareMapper());
	   }    

	   fillCrystals(diskOriginToCrystalOrigin); //See note in DiskCalorimeterMaker
      }
      


     
      // take the crystals from the CrystalMap, and keep only those who are in the annulus
      void Disk::fillCrystals(const CLHEP::Hep3Vector &crystalOriginInDisk)
      {   
	  int nRingsMax = int(2*_radiusOut/_cellSize);

          int nCrystal(0);
 	  for (int i=0;i<_crystalMap->nCrystalMax(nRingsMax);++i)
	  {
	      CLHEP::Hep2Vector xy  = _cellSize*_crystalMap->xyFromIndex(i);
              CLHEP::Hep3Vector posFF(xy.x(),xy.y(),0);
	      CLHEP::Hep3Vector pos = posFF + crystalOriginInDisk;
	      

 	      if ( !isInsideDisk(xy.x(),xy.y()) ) {_mapToCrystal.push_back(-1); continue;}

	      _crystalToMap.push_back(i);
	      _mapToCrystal.push_back(nCrystal);
	      _crystalList.push_back( Crystal(nCrystal,_id, pos, posFF) );		

	      ++nCrystal;
	  }
      }


      bool Disk::isInsideDisk(double x, double y) const
      {    	 
	 
          for (int i=1;i<_crystalMap->nApex();++i)
	  {  
              CLHEP::Hep2Vector p1(x + _cellSize*_crystalMap->apexX(i-1), y + _cellSize*_crystalMap->apexY(i-1));
              CLHEP::Hep2Vector p2(x + _cellSize*_crystalMap->apexX(i),   y + _cellSize*_crystalMap->apexY(i));

      	      //check distances. Note that the farthest distance is always at an apex in our case 
              if (calcDistToSide(p1,p2) < _radiusIn)        return false;
              if (std::max(p1.mag(),p2.mag()) > _radiusOut) return false;      
          }

          return true;
      }

      double Disk::calcDistToSide(const CLHEP::Hep2Vector &P1, const CLHEP::Hep2Vector &P0) const
      {	  
	  CLHEP::Hep2Vector v = P1-P0;
       
	  double t = -1.0*P0*v/(v*v);
	  if ( t < 0.0 ) return P0.mag();
	  if ( t > 1.0 ) return P1.mag();
	  return (P0+t*v).mag();
      }




      int Disk::idxFromPosition(double x, double y) const 
      {
           unsigned int mapIdx = _crystalMap->indexFromXY(x/_cellSize,y/_cellSize);
           if (mapIdx > _mapToCrystal.size() ) return -1;
	   return _mapToCrystal.at(mapIdx);
      }





      //find the local indexes of the crystal neighbors for a given level (level = number of rings away)
      std::vector<int> Disk::findLocalNeighbors(int crystalId, int level) const
      {
	   std::vector<int> list; 
	   std::vector<int> temp(_crystalMap->neighbors(_crystalToMap.at(crystalId),level));

	   for (unsigned int i=0;i<temp.size();++i)
	     if (_mapToCrystal.at(temp[i])>-1) list.push_back(_mapToCrystal.at(temp[i]));
	   
	   return list;
      }

      
      
      //Slightly inefficient but robust integration. Divide the area between the disk and 
      //the first few crystals into tiny squares, and sum square surface in the empty space
      //Use symmetry, do it for a quarter slice
      double Disk::estimateEmptySpace(void) const
      {
	    double sum(0),delta(0.1);
	    double RadiusMax = _radiusIn + 0.2*(_radiusOut-_radiusIn);
	    
	    for (double x=delta/2.0; x <= RadiusMax; x+=delta)
	    {
               double y0 = (x < _radiusIn) ? sqrt(_radiusIn*_radiusIn - x*x) + delta/2.0  : delta/2.0;
	       double ymax = sqrt(RadiusMax*RadiusMax-x*x);

	       for (double y=y0; y <= ymax; y+=delta)
	       {
		  int mapIdx = _crystalMap->indexFromXY( x/_cellSize,y/_cellSize);
		  int iCry   = _mapToCrystal.at(mapIdx);
		  if (iCry==-1) sum+=delta*delta;		 
	       }  
	    }

	    return 4.0*sum;
      }


      void Disk::print(std::ostream &os) const 
      {        
	CaloSection::print();
	os<<std::endl;
	os<<"Number of crystals = "<<_crystalList.size()<<std::endl;
	os<<"Radius In / Out "<<_radiusIn<<" / "<<_radiusOut<<std::endl;
      }

 
    
}

