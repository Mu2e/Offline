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
      geomUtil_(disks_,caloInfo_,geomInfo_,fullCrystalList_)
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
    
    int DiskCalorimeter::crystalIdxFromPosition(const CLHEP::Hep3Vector &pos) const 
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
    
 
    void DiskCalorimeter::print(std::ostream &os) const 
    {
       os<<"Disk calorimeter "<<std::endl;
       os<<"Number of disks :"<< disks_.size()<<std::endl;
       for (unsigned int idisk=0;idisk<disks_.size();++idisk) disk(idisk).print(os);
    }


}
