#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/CaloGeomUtil.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
#include <algorithm>


namespace mu2e {


    DiskCalorimeter::DiskCalorimeter() : 
      disks_(),
      caloInfo_(),
      geomInfo_(),
      fullCrystalList_(),  
      geomUtil_(disks_, caloInfo_, geomInfo_, fullCrystalList_)
    {}



    std::vector<int> DiskCalorimeter::neighborsByLevel(int crystalId, int level, bool rawMap)  const
    {
        int iv = fullCrystalList_.at(crystalId)->diskId();

	int offset(0);
	for (int i=0;i<iv;++i) offset += disk(i).nCrystals();

        std::vector<int> list = disk(iv).findLocalNeighbors(fullCrystalList_.at(crystalId)->localId(),level,rawMap);
	if (rawMap) transform(list.begin(), list.end(), list.begin(), [=](int i){return i>-1? i+offset :-1;});
       // else     transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  
        else     transform(list.begin(), list.end(), list.begin(),[=](int i){return i+offset;});  

	return list;
    }
    
    
    int DiskCalorimeter::crystalIdxFromPosition(const CLHEP::Hep3Vector& pos) const 
    {   
	int offset(0);
	for (unsigned int idisk=0;idisk<disks_.size();++idisk)
	{
 	   if ( geomUtil_.isInsideSection(idisk,pos) )
	   {
		 CLHEP::Hep3Vector posInSection = geomUtil_.mu2eToDisk(idisk, pos);
		 return offset + disk(idisk).idxFromPosition(posInSection.x(),posInSection.y());
	   } 
	   offset += disk(idisk).nCrystals();
	}       	  
	return -1;
    }
    

    int DiskCalorimeter::nearestIdxFromPosition(const CLHEP::Hep3Vector& pos) const 
    {                   
        int offset(0);
        CLHEP::Hep3Vector posInSection = geomUtil_.mu2eToDisk(0,pos);
        
        
        auto minDistIter = std::min_element(disks_.begin(),disks_.end(), [&](const auto& disk1, const auto& disk2)
                           { return this->deltaZ(disk1->geomInfo().origin(),pos) < this->deltaZ(disk2->geomInfo().origin(),pos); });

        for (int i=0;i<std::distance(disks_.begin(),minDistIter);++i) offset += disk(i).nCrystals();
                
        std::vector<int> cand = disks_[0]->nearestIdxFromPosition(posInSection.x(),posInSection.y());
        auto bestCandIter = std::min_element(cand.begin(),cand.end(),[&](int ic1, int ic2) {return this->deltaPerp(ic1,pos) < this->deltaPerp(ic2,pos);});

	return *bestCandIter+offset;
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
       for (unsigned int idisk=0;idisk<disks_.size();++idisk) disk(idisk).print(os);
    }


}
