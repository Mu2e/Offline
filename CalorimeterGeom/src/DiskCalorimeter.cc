//
// Geometry and identifier info about the Calorimeter.
//
// Original author B. Echenard
//

// C++ includes
#include <iostream>
#include <algorithm>

// Mu2e includes
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"

//other includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {


   bool DiskCalorimeter::isInsideCalorimeter(const CLHEP::Hep3Vector &pos) const 
   {   
       for (unsigned int idisk=0;idisk<_nSections;++idisk) if (isInsideSection(idisk, pos)) return true;	
       return false;    
   }

   bool DiskCalorimeter::isInsideSection(int idisk, const CLHEP::Hep3Vector &pos) const 
   {   
	CLHEP::Hep3Vector posInSection = toSectionFrameFF(idisk, pos);
	double posZ = posInSection.z();
	double zlim = 2*_caloGeomInfo.crystalHalfLength()+1e-6;


	if ( posZ < -1e-6 || posZ > zlim )                                          return false;
	if ( disk(idisk).idxFromPosition(posInSection.x(),posInSection.y()) == -1)  return false;	

	return true;
    }
        
    
    bool DiskCalorimeter::isContainedSection(const CLHEP::Hep3Vector &front, const CLHEP::Hep3Vector &back) const 
    {   
	//for (unsigned int idisk=0;idisk<_nSections;++idisk) 
	//   if (isInsideSection(idisk,start) && isInsideSection(idisk,stop)) return true;
	//return false;

	//this is a more efficient way of doing this in the case of two disks
	for (unsigned int idisk=0;idisk<_nSections;++idisk) 
	{
	  if ( abs(front.z()-back.z()) < (*_sections.at(idisk)).size().z() ) return true;
	}
	return false;   
    }
	 

    std::vector<int> DiskCalorimeter::neighborsByLevel(int crystalId, int level)  const
    {
	int iv = _fullCrystalList.at(crystalId)->sectionId();

	int offset(0);
	for (int i=0;i<iv;++i) offset += disk(i).nCrystals();

	std::vector<int> list = disk(iv).findLocalNeighbors( _fullCrystalList.at(crystalId)->localId() ,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  

	return list;
    }
    
    
    
    
    int DiskCalorimeter::crystalIdxFromPosition(const CLHEP::Hep3Vector &pos) const 
    {   
	int offset(0);
	for (unsigned int idisk=0;idisk<_nSections;++idisk)
	{
 	   if ( isInsideSection(idisk,pos) )
	   {
		 CLHEP::Hep3Vector posInSection = toSectionFrame(idisk, pos);
		 return offset + disk(idisk).idxFromPosition(posInSection.x(),posInSection.y());
	   } 
	   offset += disk(idisk).nCrystals();
	}       	  
	return -1;
    }
    
 
    void DiskCalorimeter::print(std::ostream &os) const 
    {
       os<<"Disk calorimeter "<<std::endl;
       os<<"Number of disks :"<< _nSections<<std::endl;
       for (unsigned int idisk=0;idisk<_nSections;++idisk) disk(idisk).print(os);
    }


}
