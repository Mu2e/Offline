//
// Create a disk and fills it with crystals.
// We assume that the real crystal position is close to the ideal one (at most a few mm differences).
// If not, the crystal navigation needs to be rewritten.
//
// Original author B Echenard
//

#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/SquareMapper.hh"
#include "Offline/CalorimeterGeom/inc/SquareShiftMapper.hh"
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <algorithm>
#include <map>

namespace mu2e {


      Disk::Disk(int id, double rCrystalIn, double rCrystalOut, double nominalCellSize,
                 double nominalCellLength, int offset, const CLHEP::Hep3Vector& diskOriginToCrystalOrigin) :
        id_(id),
        crystalList_(),
        geomInfo_(),
        radiusInCrystal_(rCrystalIn),
        radiusOutCrystal_(rCrystalOut),
        nominalCellSize_(nominalCellSize),
        globalCrystalOffset_(offset),
        mapToCrystal_(),
        crystalToMap_(),
        rowMax_(0)
      {
          crystalMap_ = std::shared_ptr<CrystalMapper>(new SquareShiftMapper());

          fillCrystalsIdeal(diskOriginToCrystalOrigin, nominalCellLength); //See note in DiskCalorimeterMaker
          //fillCrystals(diskOriginToCrystalOrigin);    //See note in DiskCalorimeterMaker
      }



      //-----------------------------------------------------------------------------
      // take the crystals from the ideal map, keeping only the crystals inside the annulus
      void Disk::fillCrystalsIdeal(const CLHEP::Hep3Vector& crystalOriginInDisk, double nominalCellLength)
      {
         int nCrystal(0);
         int nRingsMax   = int(1.5*radiusOutCrystal_/nominalCellSize_);
         int nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);

         mapToCrystal_.insert(mapToCrystal_.begin(), nCrystalMap, invalidID_);

         for (int mapIdx=0;mapIdx<nCrystalMap;++mapIdx)
         {
             CLHEP::Hep2Vector xy  = nominalCellSize_*crystalMap_->xyFromIndex(mapIdx);
             if (!isInsideDisk(xy.x(),xy.y(),nominalCellSize_,nominalCellSize_))  continue;

             //these crystals have been manually removed from the map by hand....
             if (std::abs(xy.x()-257.25) < 1.0  &&  std::abs(xy.y()-583.1) < 1.0) continue;
             if (std::abs(xy.x()-257.25) < 1.0  &&  std::abs(xy.y()+583.1) < 1.0) continue;
             if (std::abs(xy.x()+257.25) < 1.0  &&  std::abs(xy.y()-583.1) < 1.0) continue;
             if (std::abs(xy.x()+257.25) < 1.0  &&  std::abs(xy.y()+583.1) < 1.0) continue;

             CLHEP::Hep3Vector size(nominalCellSize_,nominalCellSize_,nominalCellLength);
             CLHEP::Hep3Vector posFF(xy.x(),xy.y(),0);
             CLHEP::Hep3Vector pos = posFF + crystalOriginInDisk;

             rowMax_ = std::max(rowMax_,std::abs(crystalMap_->rowFromIndex(mapIdx)));

             mapToCrystal_[mapIdx] = nCrystal;
             crystalToMap_.push_back(mapIdx);
             crystalList_.push_back(Crystal(nCrystal, id_, pos, pos, size));

             ++nCrystal;
         }
      }

      //-----------------------------------------------------------------------------
      // fill the crystal from the real position with the measured dimensions
      void Disk::fillCrystals(const CLHEP::Hep3Vector& crystalOriginInDisk)
      {
         /*
         int nRingsMax   = int(1.5*radiusOutCrystal_/nominalCellSize_);
         int nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);
         mapToCrystal_.insert(mapToCrystal_.begin(), nCrystalMap, invalidID_);

         //yes, this will be correctly done later
         WhateverDatabaseReader reader; // <-- this is the link to the database

         for (int i=0;i<reader.getNCrystals();++i)
         {
             // get position and width here
             CLHEP::Hep3Vector posFF = reader.position(i);
             CLHEP::Hep3Vector pos   = posFF + crystalOriginInDisk;
             CLHEP::Hep3Vector size  = reader.size(i);

             int mapIdx = crystalMap_->indexFromXY(pos.x()/nominalCellSize_,pos.y()/nominalCellSize_);
             CLHEP::Hep2Vector idealPosition = nominalCellSize_*crystalMap_->xyFromIndex(mapIdx);

             mapToCrystal_[mapIdx] = i;
             crystalToMap_.push_back(mapIdx);
             crystalList_.push_back( Crystal(nCrystal, id_, pos, idealPosition, size) );
         }

         fixCrystalPosition();
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

             // shortest distance between the segment P0-P1 and the origin.
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
      //check crystal overlaps (row after row) and move crystal around to fix boundary crossings
      void Disk::fixCrystalPosition()
      {
         std::map<int,std::vector<size_t>> rowToCrystalId;
         for (size_t cryId=0;cryId<crystalList_.size();++cryId){
            int irow = crystalMap_->rowFromIndex(crystalToMap_[cryId]);
            rowToCrystalId[irow].push_back(cryId);
         }

         double ymaxPrevious(-1e6);
         for (auto& kv : rowToCrystalId){
           auto& cryList = kv.second;

           //order the crystal ids by increasing x coordinate in a given row
           auto sortFunctor = [this](const auto i, const auto j)
                                  {return crystalList_[i].localPosition().x() <
                                            crystalList_[j].localPosition().x();};
           std::sort(std::begin(cryList),std::end(cryList),sortFunctor);

           //check if the crystals overlap in Y with the row below and fix the position
           double ymax(-1e6);
           for (const auto& cryIdx : cryList) {
              double crYmin    = crystalList_[cryIdx].localPosition().y()-crystalList_[cryIdx].size().y()/2.0;
              double tolerance = crYmin - ymaxPrevious;
              if (tolerance < 0) {
                 auto newPosition = crystalList_[cryIdx].localPosition() - CLHEP::Hep3Vector(0,tolerance,0);
                 crystalList_[cryIdx].setLocalPosition(newPosition);
              }
              ymax = std::max(ymax,crystalList_[cryIdx].localPosition().y()+crystalList_[cryIdx].size().y()/2.0);
           }
           ymaxPrevious = ymax;

           //check if the crystal overlap in X, and move them to the right if needed
           //double origXend = crystalList_[cryList.back()].localPosition().x();
           for (size_t i=1;i<cryList.size();++i) {
              size_t idx0 = cryList[i-1];
              size_t idx1 = cryList[i];
              double crXmax0 = crystalList_[idx0].localPosition().x()+crystalList_[idx0].size().x()/2.0;
              double crXmin1 = crystalList_[idx1].localPosition().x()-crystalList_[idx1].size().x()/2.0;
              double tolerance = crXmin1 - crXmax0;
              if (tolerance < 0) {
                 auto newPosition = crystalList_[idx1].localPosition() - CLHEP::Hep3Vector(tolerance,0,0);
                 crystalList_[idx1].setLocalPosition(newPosition);
              }
           }

           //recenter all crystals by shifting them to the left if we had to adjust some
           //double newXend = crystalList_[cryList.back()].localPosition().x();
           //double shiftPerCrystal = (newXend-origXend)/2.0;
           //if (shiftPerCrystal>0){
           //  for (const auto& cryIdx : cryList) {
           //      auto newPosition = crystalList_[cryIdx].localPosition() - CLHEP::Hep3Vector(shiftPerCrystal,0,0);
           //      crystalList_[cryIdx].setLocalPosition(newPosition);
           //   }
           //}

           // check that all crystals are still in the disk enveloppe
           for (const auto& cryIdx : cryList) {
             auto position = crystalList_[cryIdx].localPosition();
             auto size     = crystalList_[cryIdx].size();
             if (!isInsideDisk(position.x(),position.y(),size.x(),size.y()))
               throw cet::exception("Disk") << " The crystal at position="<<position<<" does not fit inside the disk...\n";
           }
         }
      }







      //-----------------------------------------------------------------------------
      //calculate the bounding box surrounding a row of crystals
      void Disk::boundingBoxes(int thisRow, std::vector<double>& params) const
      {
         params.clear();

         std::vector<size_t> cryList;
         for (size_t i=0;i<crystalList_.size();++i){
           int irow = crystalMap_->rowFromIndex(crystalToMap_[i]);
           if (irow == thisRow) cryList.push_back(i);
         }
         if (cryList.empty()) return;

         auto sortFunctor = [this](const auto i, const auto j)
                                  {return crystalList_[i].localPosition().x() <
                                          crystalList_[j].localPosition().x();};
         std::sort(std::begin(cryList),std::end(cryList),sortFunctor);

         double ymin(1e6),ymax(-1e6),xInmin(-1e6),xImax(1e6);
         for (size_t i=0;i<cryList.size();++i) {
           ymin = std::min(ymin,crystalList_[cryList[i]].localPosition().y()-crystalList_[cryList[i]].size().y()/2.0);
           ymax = std::max(ymax,crystalList_[cryList[i]].localPosition().y()+crystalList_[cryList[i]].size().y()/2.0);
           if (i==0) continue;
           double x0 = crystalList_[cryList[i-1]].localPosition().x()+crystalList_[cryList[i-1]].size().x()/2.0;
           double x1 = crystalList_[cryList[i]].localPosition().x()  -crystalList_[cryList[i]].size().x()/2.0;
           if (x1-x0 > 2*crystalList_[cryList[i]].size().x()){xInmin = x0; xImax=x1;}
         }

         double xOutmin = crystalList_[cryList.front()].localPosition().x()-crystalList_[cryList.front()].size().x()/2.0;
         double xOutmax = crystalList_[cryList.back()].localPosition().x() +crystalList_[cryList.back()].size().x()/2.0;
         params.push_back(xOutmin);
         params.push_back(xOutmax);
         params.push_back(ymin);
         params.push_back(ymax);
         if (xInmin>-1e5) {params.push_back(xInmin);params.push_back(xImax);}
      }


      //-----------------------------------------------------------------------------
      int Disk::idMinCrystalInside(int row) const
      {
         if (std::abs(row) > rowMax_) return invalidID_;
         int idx(0);
         while (mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)] == invalidID_) ++idx;
         return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];
      }

      //-----------------------------------------------------------------------------
      int Disk::idMaxCrystalInside(int row) const
      {
         if (std::abs(row) > rowMax_) return invalidID_;
         int idx = int(radiusOutCrystal_/nominalCellSize_)+2;
         while ( !isMapIdxValid(crystalMap_->indexFromRowCol(row,idx)) ) --idx;
         while ( mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)]==invalidID_ && idx>0) --idx;
         return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];
      }



      //-----------------------------------------------------------------------------
      int Disk::idxFromPosition(double x, double y) const
      {
         // this only works for ideal crystal placement
         //int mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);
         //if (mapIdx < mapToCrystal_.size() ) return mapToCrystal_.at(mapIdx);
         //return invalidID_;

         int mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);

         if (isCrystalIdxValid(mapIdx) && isInsideCrystal(mapToCrystal_[mapIdx],x,y)) return mapToCrystal_[mapIdx];

         const int level(1);
         const auto neighbors(crystalMap_->neighbors(mapIdx,level));
         for (const auto& idx : neighbors) {
            if (isCrystalIdxValid(idx) && isInsideCrystal(mapToCrystal_[idx],x,y)) return mapToCrystal_[idx];
         }

         return invalidID_;
      }




      //-----------------------------------------------------------------------------
      bool Disk::isInsideCrystal(int icry, double x, double y) const
      {
         if (icry >= static_cast<int>(crystalList_.size())) return false;
         return crystalMap_->isInsideCrystal(x, y, crystalList_[icry].localPosition(), crystalList_[icry].size());
      }


      //-----------------------------------------------------------------------------
      //find the local indexes of the crystal neighbors for a given level (level = number of rings away)
      std::vector<int> Disk::findLocalNeighbors(int crystalId, int level) const
      {
         std::vector<int> list;
         std::vector<int> temp(crystalMap_->neighbors(crystalToMap_.at(crystalId),level));

         for (const auto& mapIdx : temp) {
           if (isCrystalIdxValid(mapIdx)) list.push_back(mapToCrystal_.at(mapIdx));
         }
         return list;
      }


      //-----------------------------------------------------------------------------
      //find the nearest crystals from the position based on ideal mapping
      std::vector<int> Disk::nearestNeighborsFromPos(double x, double y) const
      {
         int level(1);
         std::vector<int> list;

         auto mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);
         if (isCrystalIdxValid(mapIdx)) list.push_back(mapToCrystal_.at(mapIdx));

         while (list.empty()) {
           for (auto mapIdx : crystalMap_->neighbors(mapIdx,level)) {
             if (isCrystalIdxValid(mapIdx)) list.push_back(mapToCrystal_[mapIdx]);
           }
           ++level;
         }
         return list;
      }



     //-----------------------------------------------------------------------------
     const bool Disk::isCrystalIdxValid(int i) const {
       return i < static_cast<int>(mapToCrystal_.size()) && mapToCrystal_[i]!=invalidID_;
     }

     //-----------------------------------------------------------------------------
     const bool Disk::isMapIdxValid(int i) const {
       return i < static_cast<int>(mapToCrystal_.size());
     }


     //-----------------------------------------------------------------------------
     void Disk::print(std::ostream &os) const
     {
         os<<"Disk                  "<<id_<<std::endl;
         os<<"Number of crystals    "<<crystalList_.size()<<std::endl;
         os<<"Radius In             "<<geomInfo_.innerEnvelopeR()<<std::endl;
         os<<"Radius Out            "<<geomInfo_.outerEnvelopeR()<<std::endl;
         os<<"origin                "<<geomInfo_.originLocal()<<std::endl;
         os<<"origin Mu2e           "<<geomInfo_.origin()<<std::endl;
         os<<"size                  "<<geomInfo_.size()<<std::endl;
         os<<"rotation              "<<geomInfo_.rotation()<<std::endl;
         os<<"originToCrystalOrigin "<<geomInfo_.originToCrystalOrigin()<<std::endl;
         os<<"z Front               "<<geomInfo_.frontFaceCenter().z()<<std::endl;
         os<<"z Back                "<<geomInfo_.backFaceCenter().z()<<std::endl;
         os<<"r In tracker          "<<geomInfo_.innerEnvelopeR()<<std::endl;
         os<<"r Out tracker         "<<geomInfo_.outerEnvelopeR()<<std::endl;
     }

}

