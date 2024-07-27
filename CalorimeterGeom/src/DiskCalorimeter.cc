#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
#include <algorithm>


namespace mu2e {


    DiskCalorimeter::DiskCalorimeter() :
      disks_(),
      fullCrystalList_(),
      caloInfo_(),
      geomUtil_(disks_, fullCrystalList_)
    {}


    const std::vector<int>& DiskCalorimeter::neighbors(int crystalId) const
    {
        return fullCrystalList_.at(crystalId)->neighbors();
    }

    const std::vector<int>& DiskCalorimeter::nextNeighbors(int crystalId) const
    {
        return fullCrystalList_.at(crystalId)->nextNeighbors();
    }

    std::vector<int> DiskCalorimeter::neighborsByLevel(int crystalId, int level)  const
    {
        int iv = fullCrystalList_.at(crystalId)->diskID();
        int offset = disk(iv).crystalOffset();

        std::vector<int> list = disk(iv).findLocalNeighbors(fullCrystalList_.at(crystalId)->localID(),level);
        transform(list.begin(), list.end(), list.begin(),[=](int i){return (i<0) ? i : i+offset;});

        return list;
    }


    int DiskCalorimeter::crystalIdxFromPosition(const CLHEP::Hep3Vector& pos) const
    {
        for (size_t idisk=0;idisk<disks_.size();++idisk)
        {
           if ( geomUtil_.isInsideSection(idisk,pos) )
           {
                 CLHEP::Hep3Vector posInSection = geomUtil_.mu2eToDisk(idisk, pos);
                 return disk(idisk).crystalOffset() + disk(idisk).idxFromPosition(posInSection.x(),posInSection.y());
           }
        }
        return -1;
    }


    int DiskCalorimeter::nearestIdxFromPosition(const CLHEP::Hep3Vector& pos) const
    {
        CLHEP::Hep3Vector posInSection = geomUtil_.mu2eToDisk(0,pos);

        auto minPred     = [&](const auto& disk1, const auto& disk2)
                           {return this->deltaZ(disk1->geomInfo().origin(),pos) < this->deltaZ(disk2->geomInfo().origin(),pos);};
        auto minDistIter = std::min_element(disks_.begin(),disks_.end(),minPred);

        std::vector<int> cand = disks_[0]->nearestNeighborsFromPos(posInSection.x(),posInSection.y());
        auto bestPred     = [&](int ic1, int ic2) {return this->deltaPerp(ic1,pos) < this->deltaPerp(ic2,pos);};
        auto bestCandIter = std::min_element(cand.begin(),cand.end(),bestPred);

        return *bestCandIter+(*minDistIter)->crystalOffset();
    }

    double DiskCalorimeter::deltaZ(const CLHEP::Hep3Vector& p1, const CLHEP::Hep3Vector& p2) const
    {
        return std::abs(p1.z()-p2.z());
    }

    double DiskCalorimeter::deltaPerp(int ic, const CLHEP::Hep3Vector& pos) const
    {
        return sqrt((fullCrystalList_.at(ic)->position()-pos).perp2());
    }



    void DiskCalorimeter::print(std::ostream &os) const
    {
       os<<"Disk calorimeter "<<std::endl;
       os<<"Number of disks :"<< disks_.size()<<std::endl;
       for (size_t idisk=0;idisk<disks_.size();++idisk) disk(idisk).print(os);
    }


}
