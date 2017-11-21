//
// Create a disk and fills it with crystals
//
// Original author B Echenard
//

// Notes: CrystalMap tesselates a plane with square crystals and let you know which 
//        ones are fully contained inside an annulus

#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/SquareMapper.hh"
#include "CalorimeterGeom/inc/SquareShiftMapper.hh"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

      
      Disk::Disk(int id, double rin, double rout,double cellSize, bool shiftCrystal,
		 const CLHEP::Hep3Vector& diskOriginToCrystalOrigin) : 
	crystalList_(), 
	id_(id), 
        geomInfo_(),
	radiusIn_(rin),
	radiusOut_(rout),
	cellSize_(cellSize),
	mapToCrystal_(),
	crystalToMap_()
      { 
	   geomInfo_.originToCrystalOrigin(diskOriginToCrystalOrigin);
           
           if (shiftCrystal) crystalMap_ = std::shared_ptr<CrystalMapper>(new SquareShiftMapper());
	   else              crystalMap_ = std::shared_ptr<CrystalMapper>(new SquareMapper());

	   fillCrystals(diskOriginToCrystalOrigin); //See note in DiskCalorimeterMaker
      }
      

     
      //-----------------------------------------------------------------------------
      // take the crystals from the CrystalMap, and keep only those who are in the annulus
      void Disk::fillCrystals(const CLHEP::Hep3Vector& crystalOriginInDisk)
      {   
	  int nRingsMax = int(2*radiusOut_/cellSize_);

          int nCrystal(0);
 	  for (int i=0;i<crystalMap_->nCrystalMax(nRingsMax);++i)
	  {
	      CLHEP::Hep2Vector xy  = cellSize_*crystalMap_->xyFromIndex(i);	      
 	      if ( !isInsideDisk(xy.x(),xy.y()) ) {mapToCrystal_.push_back(-1); continue;}
 	      
              //these crystals have been manually removed from the map for whatever stupid reason...
              if ( std::abs(xy.x()-257.25) < 1.0 && std::abs(xy.y()-583.1) < 1.0) {mapToCrystal_.push_back(-1); continue;}
              if ( std::abs(xy.x()-257.25) < 1.0 && std::abs(xy.y()+583.1) < 1.0) {mapToCrystal_.push_back(-1); continue;}
              if ( std::abs(xy.x()+257.25) < 1.0 && std::abs(xy.y()-583.1) < 1.0) {mapToCrystal_.push_back(-1); continue;}
              if ( std::abs(xy.x()+257.25) < 1.0 && std::abs(xy.y()+583.1) < 1.0) {mapToCrystal_.push_back(-1); continue;}

              CLHEP::Hep3Vector posFF(xy.x(),xy.y(),0);
	      CLHEP::Hep3Vector pos = posFF + crystalOriginInDisk;
              
              crystalToMap_.push_back(i);
	      mapToCrystal_.push_back(nCrystal);
	      crystalList_.push_back( Crystal(nCrystal,id_, pos) );		
	      ++nCrystal;
	  }            
      }


      //-----------------------------------------------------------------------------
      int Disk::idMinCrystalInside(int row)
      {
          int idx(0);
          while (mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)]<0) ++idx;
         
          return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];      
      }

      int Disk::idMaxCrystalInside(int row)
      {
          int idx = int(radiusOut_/cellSize_)+2;
          while (mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)]<0 && idx>0) --idx;
          return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];      
      }

      //-----------------------------------------------------------------------------
      bool Disk::isInsideDisk(double x, double y) const
      {    	 	 
          for (int i=1;i<crystalMap_->nApex();++i)
	  {  
              CLHEP::Hep2Vector p1(x + cellSize_*crystalMap_->apexX(i-1), y + cellSize_*crystalMap_->apexY(i-1));
              CLHEP::Hep2Vector p2(x + cellSize_*crystalMap_->apexX(i),   y + cellSize_*crystalMap_->apexY(i));

      	      //check distances. Note that the farthest distance is always at an apex in our case 
              if (calcDistToSide(p1,p2) < radiusIn_)        return false;
              if (std::max(p1.mag(),p2.mag()) > radiusOut_) return false;      
          }

          return true;
      }

      //-----------------------------------------------------------------------------
      double Disk::calcDistToSide(const CLHEP::Hep2Vector &P1, const CLHEP::Hep2Vector &P0) const
      {	  
	  CLHEP::Hep2Vector v = P1-P0;
       
	  double t = -1.0*P0*v/(v*v);
	  if ( t < 0.0 ) return P0.mag();
	  if ( t > 1.0 ) return P1.mag();
	  return (P0+t*v).mag();
      }




     //-----------------------------------------------------------------------------
     int Disk::idxFromPosition(double x, double y) const 
      {
           unsigned int mapIdx = crystalMap_->indexFromXY(x/cellSize_,y/cellSize_);
           if (mapIdx > mapToCrystal_.size() ) return -1;
	   return mapToCrystal_.at(mapIdx);
      }





     //-----------------------------------------------------------------------------
     //find the local indexes of the crystal neighbors for a given level (level = number of rings away)
     std::vector<int> Disk::findLocalNeighbors(int crystalId, int level, bool raw) const
     {
          std::vector<int> list; 
	  std::vector<int> temp(crystalMap_->neighbors(crystalToMap_.at(crystalId),level));

          for (unsigned int i=0;i<temp.size();++i)
          {
	     if (raw) {list.push_back(mapToCrystal_.at(temp[i]));}
             else {if (mapToCrystal_.at(temp[i])>-1) list.push_back(mapToCrystal_.at(temp[i]));}	      
          } 

          return list;
     }      


     //-----------------------------------------------------------------------------
     //find the nearest crystals from the position
     std::vector<int> Disk::nearestIdxFromPosition(double x, double y) const 
     {
          int level(1);
          std::vector<int> list;

          unsigned int mapIdx = crystalMap_->indexFromXY(x/cellSize_,y/cellSize_);
          if (mapIdx < mapToCrystal_.size() && mapToCrystal_.at(mapIdx)>-1) list.push_back( mapToCrystal_.at(mapIdx) );

          while(list.empty())
          {
             std::vector<int> temp(crystalMap_->neighbors(mapIdx,level));
             for (auto it : temp) 
                if (mapToCrystal_.at(it)>-1) list.push_back(mapToCrystal_.at(it));
             ++level;                  
          }

	  return list;
     }


     //Slightly inefficient but robust integration. Divide the area between the disk and 
     //the first few crystals into tiny squares, and sum square surface in the empty space
     //Use symmetry, do it for a quarter slice
     double Disk::estimateEmptySpace(void) const
     {
	   double sum(0),delta(0.1);
	   double RadiusMax = radiusIn_ + 0.2*(radiusOut_-radiusIn_);

	   for (double x=delta/2.0; x <= RadiusMax; x+=delta)
	   {
              double y0 = (x < radiusIn_) ? sqrt(radiusIn_*radiusIn_ - x*x) + delta/2.0  : delta/2.0;
	      double ymax = sqrt(RadiusMax*RadiusMax-x*x);

	      for (double y=y0; y <= ymax; y+=delta)
	      {
		 int mapIdx = crystalMap_->indexFromXY( x/cellSize_,y/cellSize_);
		 int iCry   = mapToCrystal_.at(mapIdx);
		 if (iCry==-1) sum+=delta*delta;		 
	      }  
	   }

	   return 4.0*sum;
     }


     void Disk::print(std::ostream &os) const
     {
          os<<"Disk                   "<<id_<<std::endl;
    	  os<<"Number of crystals =   "<<crystalList_.size()<<std::endl;
	  os<<"Radius In / Out        "<<radiusIn_<<" / "<<radiusOut_<<std::endl;
          os<<"origin                 "<<geomInfo_.originLocal()<<std::endl;
          os<<"origin Mu2e            "<<geomInfo_.origin()<<std::endl;
          os<<"size                   "<<geomInfo_.size()<<std::endl;
          os<<"rotation               "<<geomInfo_.rotation()<<std::endl;
          os<<"originToCrystalOrigin  "<<geomInfo_.originToCrystalOrigin()<<std::endl;
          os<<"z Front                "<<geomInfo_.frontFaceCenter().z()<<std::endl;
          os<<"z Back                 "<<geomInfo_.backFaceCenter().z()<<std::endl;
          os<<"r In tracker           "<<geomInfo_.innerEnvelopeR()<<std::endl;
          os<<"r Out tracker          "<<geomInfo_.outerEnvelopeR()<<std::endl;
     }  
    
}

