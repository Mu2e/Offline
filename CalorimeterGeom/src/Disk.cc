//
// Create a disk and fills it with crystals.
// We assume that the real crystal position is close to the ideal one (at most a few mm differences).
// If not, the crystal navigation needs to be rewritten
//
// Original author B Echenard
//
// Notes: CrystalMap tesselates a plane with square crystals and let you know which
//        ones are fully contained inside an annulus
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
          int nRingsMax   = int(1.5*radiusOutCrystal_/nominalCellSize_);
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
              crystalList_.push_back(Crystal(nCrystal, id_, pos, pos, size));
              ++nCrystal;
          }

          fixCrystalPosition();
      }

      //-----------------------------------------------------------------------------
      // fill the crystal from the real position with the measured dimensions
      void Disk::fillCrystals(const CLHEP::Hep3Vector& crystalOriginInDisk)
      {
          /*
          int nCrystal(0);
          int nRingsMax   = int(1.5*radiusOutCrystal_/nominalCellSize_);
          int nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);

          mapToCrystal_.insert(mapToCrystal_.begin(), nCrystalMap, -1);

          //yes, this will be correctly done later
          WhateverDatabaseReader reader; // <-- this is the link to the database

          for (int i=0;i<nCrystalMap;++i)
          {
              // get position and width here
              CLHEP::Hep3Vector posFF = reader.position(i);
              CLHEP::Hep3Vector pos   = posFF + crystalOriginInDisk;
              CLHEP::Hep3Vector size  = reader.size(i);

              if (!isInsideDisk(pos.x(),pos.y(),size.x(),size.y()))
                throw cet::exception("Disk") << " The crystal at position="<<posFF<<" does not fit inside the disk...\n";

              int mapIdx = crystalMap_->indexFromXY(pos.x()/nominalCellSize_,pos.y()/nominalCellSize_);
              CLHEP::Hep2Vector idealPosition = nominalCellSize_*crystalMap_->xyFromIndex(mapIdx);

              mapToCrystal_[mapIdx] = nCrystal;
              crystalToMap_.push_back(mapIdx);
              crystalList_.push_back( Crystal(nCrystal, id_, pos, idealPosition, size) );
              ++nCrystal;
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
      void Disk::fixCrystalPosition()
      {
          std::map<int,std::vector<int>> rowToCrystalId;
          for (size_t cryId=0;cryId<crystalList_.size();++cryId){
             int irow = crystalMap_->rowFromIndex(crystalToMap_[cryId]);
             rowToCrystalId[irow].push_back(cryId);
          }

          float ymaxPrevious(-1e6);
          for (auto& kv : rowToCrystalId){
            auto& cryList = kv.second;

            //order the crystal ids by increasing x coordinate in a given row
            auto sortFunctor = [this](const int i, const int j){return crystalList_[i].localPosition().x() <crystalList_[j].localPosition().x();};
            std::sort(std::begin(cryList),std::end(cryList),sortFunctor);

            //check if the crystals overlap in Y with the row below and fix the position
            double ymax(-1e6);
            for (const auto& cryIdx : cryList) {
               float crYmin    = crystalList_[cryIdx].localPosition().y()-crystalList_[cryIdx].size().y()/2.0;
               float tolerance = crYmin - ymaxPrevious;
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
               float crXmax0 = crystalList_[idx0].localPosition().x()+crystalList_[idx0].size().x()/2.0;
               float crXmin1 = crystalList_[idx1].localPosition().x()-crystalList_[idx1].size().x()/2.0;
               float tolerance = crXmin1 - crXmax0;
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
          }
      }



      //-----------------------------------------------------------------------------
      //calculate the bounding box surrounding a row of crystals
      void Disk::boundingBoxes(int thisRow, std::vector<double>& params)
      {
        params.clear();

        std::vector<int> cryList;
        for (size_t i=0;i<crystalList_.size();++i){
           int irow = crystalMap_->rowFromIndex(crystalToMap_[i]);
           if (irow == thisRow) cryList.push_back(i);
        }
        auto sortFunctor = [this](const int i, const int j){return crystalList_[i].localPosition().x() <crystalList_[j].localPosition().x();};
        std::sort(std::begin(cryList),std::end(cryList),sortFunctor);

        if (cryList.empty()) return;

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
      int Disk::idMinCrystalInside(int row)
      {
          int idx(0);
          while (mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)]<0) ++idx;
          return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];
      }

      int Disk::idMaxCrystalInside(int row)
      {
          int idx = int(radiusOutCrystal_/nominalCellSize_)+2;
          while (crystalMap_->indexFromRowCol(row,idx) > int(mapToCrystal_.size())) --idx;
          while (mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)]<0 && idx>0) --idx;
          return mapToCrystal_[crystalMap_->indexFromRowCol(row,idx)];
      }











      //-----------------------------------------------------------------------------
      int Disk::idxFromPosition(double x, double y) const
      {
          // this only works for ideal crystal placement
           //unsigned int mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);
          //if (mapIdx < mapToCrystal_.size() ) return mapToCrystal_.at(mapIdx);
          // return -1;

          unsigned int mapIdx = crystalMap_->indexFromXY(x/nominalCellSize_,y/nominalCellSize_);
          if (mapIdx >= mapToCrystal_.size()) return -1;

          int crystalId = mapToCrystal_.at(mapIdx);
          if (isInsideCrystal(crystalId,x,y)) return crystalId;

          std::vector<int> neighbors;
          if (crystalId > -1) neighbors = crystalList_[crystalId].neighbors();
          else {
             std::vector<int> temp(crystalMap_->neighbors(mapIdx,1));
             for (const auto& val : temp) if (mapToCrystal_.at(val) >-1) neighbors.push_back(mapToCrystal_.at(val));
          }

          for (unsigned it : neighbors)
             if (it < mapToCrystal_.size() && isInsideCrystal(mapToCrystal_[it],x,y)) return mapToCrystal_[it];

          return -1;
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

