//
// Create a disk and fills it with crystals. 
// We assume that the real crystal position is close to the ideal one. If not, 
//     this module needs to be completely rewritten
//
// Original author B Echenard
//

// Notes: CrystalMap tesselates a plane with square crystals and let you know which 
//        ones are fully contained inside an annulus

#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/SquareMapper.hh"
#include "CalorimeterGeom/inc/SquareShiftMapper.hh"
#include "CalorimeterGeom/inc/CrystalCondReader.hh"
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

      
      Disk::Disk(int id, double rin, double rout, double rCrystalIn, double rCrystalOut, double nominalCellSize, 
                 int offset, const CLHEP::Hep3Vector& diskOriginToCrystalOrigin) : 
	id_(id), 
	crystalList_(), 
        geomInfo_(),
	radiusIn_(rin),
	radiusOut_(rout),
	radiusInCrystal_(rCrystalIn),
	radiusOutCrystal_(rCrystalOut),
	nominalCellSize_(nominalCellSize),
        globalCrystalOffset_(offset),
	mapToCrystal_(),
	crystalToMap_()
      { 
	   geomInfo_.originToCrystalOrigin(diskOriginToCrystalOrigin);           
           crystalMap_ = std::shared_ptr<CrystalMapper>(new SquareShiftMapper());

	   fillCrystalsIdeal(diskOriginToCrystalOrigin); //See note in DiskCalorimeterMaker
	   //fillCrystals(diskOriginToCrystalOrigin); //See note in DiskCalorimeterMaker
      }
      

     
      //-----------------------------------------------------------------------------
      // take the crystals from the ideal map, keeping only the crystals inside the annulus
      void Disk::fillCrystalsIdeal(const CLHEP::Hep3Vector& crystalOriginInDisk)
      {   
          int nCrystal(0);
	  int nRingsMax   = int(1.5*radiusOut_/nominalCellSize_);
          int nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);
          
          mapToCrystal_.insert(mapToCrystal_.begin(), nCrystalMap, -1);

 	  for (int i=0;i<nCrystalMap;++i)
	  {
              CLHEP::Hep2Vector xy  = nominalCellSize_*crystalMap_->xyFromIndex(i);	      
 	      if (!isInsideDisk(xy.x(),xy.y(),nominalCellSize_,nominalCellSize_))  continue;
 	      
              //these crystals have been manually removed from the map by hand....
              if (std::abs(xy.x()-257.25) < 1.0  &&  std::abs(xy.y()-583.1) < 1.0) continue;
              if (std::abs(xy.x()-257.25) < 1.0  &&  std::abs(xy.y()+583.1) < 1.0) continue;
              if (std::abs(xy.x()+257.25) < 1.0  &&  std::abs(xy.y()-583.1) < 1.0) continue;
              if (std::abs(xy.x()+257.25) < 1.0  &&  std::abs(xy.y()+583.1) < 1.0) continue;

              CLHEP::Hep3Vector size(nominalCellSize_,nominalCellSize_,200.0);
              CLHEP::Hep3Vector posFF(xy.x(),xy.y(),0);
	      CLHEP::Hep3Vector pos = posFF + crystalOriginInDisk;
             
	      mapToCrystal_[i] = nCrystal;
              crystalToMap_.push_back(i);
	      crystalList_.push_back( Crystal(nCrystal, id_, pos, size) );		
	      ++nCrystal;
	  }
      }

      
      //-----------------------------------------------------------------------------
      // fil the crystal from the real position with the measured dimensions
      void Disk::fillCrystals(const CLHEP::Hep3Vector& crystalOriginInDisk)
      {                       
          /*
          int nCrystal(0);
	  int nRingsMax   = int(1.5*radiusOut_/nominalCellSize_);
          int nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);
         
          mapToCrystal_.insert(mapToCrystal_.begin(), nCrystalMap, -1);
          
          //yes, this will be correctly done later
          CrystalCondReader reader("/nfs/home/echenard/crystalCoord.txt");
                    
 	  for (int i=0;i<reader.nCrystal();++i)
	  {
	      // get position and width here
              CLHEP::Hep3Vector posFF = reader.position(i);     
	      CLHEP::Hep3Vector pos   = posFF + crystalOriginInDisk;
              CLHEP::Hep3Vector size  = reader.size(i);     
 	      
              if (!isInsideDisk(pos.x(),pos.y(),size.x(),size.y())) 
                throw cet::exception("Disk") << " The crystal at position="<<posFF<<" does not fit inside the disk...\n";
 	      	      
              int mapIdx = crystalMap_->indexFromXY(pos.x()/nominalCellSize_,pos.y()/nominalCellSize_);
              
 	      mapToCrystal_[mapIdx] = nCrystal;
              crystalToMap_.push_back(mapIdx);
	      crystalList_.push_back( Crystal(nCrystal, id_, pos, size) );		
	      ++nCrystal;             
	  } 
          */ 
      }
     



      //-----------------------------------------------------------------------------
      bool Disk::isInsideDisk(double x, double y, double widthX, double widthY) const
      {    	 	 
          std::vector<double> apexX = crystalMap_->apexX();
          std::vector<double> apexY = crystalMap_->apexY();
          
          for (size_t i=1;i<apexX.size();++i)
	  {  
              CLHEP::Hep2Vector p0(x + widthX*apexX[i-1], y + widthY*apexY[i-1]);
              CLHEP::Hep2Vector p1(x + widthX*apexX[i],   y + widthY*apexY[i]);
              
              //shortest distance between the segment P0-P1 and the origin. 
              CLHEP::Hep2Vector v = p1-p0;
	      double t = -p0*v/(v*v);              
              double mindist = (p0+t*v).mag();
              if (t < 0.0) mindist = p0.mag();
	      if (t > 1.0) mindist = p1.mag();
	      
              //farthest distance is always at an apex
              double maxdist = std::max(p0.mag(),p1.mag());
              
              if (mindist < radiusInCrystal_ || maxdist > radiusOutCrystal_) return false;
          }

          return true;
      }

      //-----------------------------------------------------------------------------
      //check that the crystals fit. If not, reduce their size to make the mapping free of overlaps
      void Disk::checkCrystalSize()
      {
          for (size_t i=0;i<crystalList_.size();++i)
          {
              Crystal center = crystalList_[i];
              auto centerSize = center.size();
              
              auto neighbors = findLocalNeighbors(i,1);
              for (size_t j=0;j<neighbors.size();++j)
              {
                 Crystal adjacent = crystalList_[j];
                 double dx = center.position().x()-adjacent.position().x()-center.size().x()-adjacent.size().x();
                 double dy = center.position().y()-adjacent.position().y()-center.size().y()-adjacent.size().y();
                 
                 //adjust crystal size if x or y is too short
                 if (dx<-1e-3 && dy <0.5*adjacent.size().y()) std::cout<<"Reduce crystal "<<i<<" size  by dx="<<dx<<std::endl;
                 if (dx<-1e-3 && dy <0.5*adjacent.size().y()) center.adjustSize(centerSize-CLHEP::Hep3Vector(dx,0,0));
                 if (dx<-1e-3 && dy <0.5*adjacent.size().y()) std::cout<<"Reduce crystal "<<i<<" size  by dy="<<dy<<std::endl;
                 if (dy<-1e-3 && dx <0.5*adjacent.size().x()) center.adjustSize(centerSize-CLHEP::Hep3Vector(0,dy,0));
              }
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
          int idx = int(radiusOut_/nominalCellSize_)+2;
          while (crystalMap_->indexFromRowCol(row,idx) > int(mapToCrystal_.size())) --idx;          
          while (mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)]<0 && idx>0) --idx;
          return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];      
      }







      


     
      //-----------------------------------------------------------------------------
      int Disk::idxFromPosition(double x, double y) const 
      {
          unsigned int mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);
          if (mapIdx < mapToCrystal_.size() ) return mapToCrystal_.at(mapIdx);
 	  return -1;

          /*
          unsigned int mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);          
          if (mapIdx >= mapToCrystal_.size()) return -1;

          int crystalId = mapToCrystal_.at(mapIdx);
          if (isInsideCrystal(crystalId,x,y)) return crystalId;

          std::vector<int> neighbors;
          if (crystalId > -1) neighbors = crystalList_[crystalId].neighbors();
          else {
             std::vector<int> temp(crystalMap_->neighbors(mapIdx,1));
             for (auto val : temp) if (mapToCrystal_.at(val) >-1) neighbors.push_back(mapToCrystal_.at(val));
          }

          for (unsigned it : neighbors)
             if (it < mapToCrystal_.size() && isInsideCrystal(mapToCrystal_[it],x,y)) return mapToCrystal_[it];

          return -1;
          */
      }

      bool Disk::isInsideCrystal(int icry, double x, double y) const
      {
          if (icry <0) return false;
          return crystalMap_->isInsideCrystal(x, y, crystalList_[icry].localPosition(), crystalList_[icry].size());
      }





      //-----------------------------------------------------------------------------
      //find the local indexes of the crystal neighbors for a given level (level = number of rings away)
      std::vector<int> Disk::findLocalNeighbors(int crystalId, int level, bool raw) const
      {
           std::vector<int> list; 
	   std::vector<int> temp(crystalMap_->neighbors(crystalToMap_.at(crystalId),level));

           for (size_t i=0;i<temp.size();++i)
           {
	      if (raw) list.push_back(mapToCrystal_.at(temp[i]));
              else 
                {if (mapToCrystal_.at(temp[i]) >-1) list.push_back(mapToCrystal_.at(temp[i]));}      
           } 

           return list;
      }      


      //-----------------------------------------------------------------------------
      //find the nearest crystals from the position based on ideal mapping
      std::vector<int> Disk::nearestIdxFromPosition(double x, double y) const 
      {
           int level(1);
           std::vector<int> list;

           unsigned mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);
           if (mapIdx < mapToCrystal_.size() && mapToCrystal_.at(mapIdx)>-1) list.push_back( mapToCrystal_.at(mapIdx) );

           while (list.size()<2)
           {
              std::vector<int> temp(crystalMap_->neighbors(mapIdx,level));
              for (unsigned it : temp) 
                 if (it < mapToCrystal_.size() && mapToCrystal_.at(it)>-1) list.push_back(mapToCrystal_.at(it));
              ++level;                  
           }

	   return list;
      }


      //-----------------------------------------------------------------------------
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
		  int mapIdx = crystalMap_->indexFromXY( x/nominalCellSize_,y/nominalCellSize_);
		  if (mapToCrystal_.at(mapIdx) <0) sum+=delta*delta;		 
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

