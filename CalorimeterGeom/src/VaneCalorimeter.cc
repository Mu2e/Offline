//
// Geometry and identifier info about the VaneCalorimeter.
//
//
// $Id: VaneCalorimeter.cc,v 1.7 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Original author R. Bernstein and Rob Kutschke
//

#include "CalorimeterGeom/inc/VaneCalorimeter.hh"

#include <algorithm>


namespace mu2e {


    bool VaneCalorimeter::isInsideCalorimeter(const CLHEP::Hep3Vector &pos) const 
    {   
        for (unsigned int ivane=0;ivane<_nSections;++ivane) if (isInsideSection(ivane,pos)) return true;
        return false;    
    }


    bool VaneCalorimeter::isInsideSection(int ivane, const CLHEP::Hep3Vector &pos) const 
    {   
	//in the Section coordinate system, the crystal are oriented along the y direction
	CLHEP::Hep3Vector posInSection = toSectionFrameFF(ivane, pos);
	double zlim = 2*_caloGeomInfo.crystalHalfLength()+1e-6;

	if (posInSection.x() < -vane(ivane).xActiveHalf() || posInSection.x() > vane(ivane).xActiveHalf()) return false;      
	if (posInSection.y() < -vane(ivane).yActiveHalf() || posInSection.y() > vane(ivane).yActiveHalf()) return false;      
	if (posInSection.z() <  -1e-6                     || posInSection.z() > zlim  )                    return false;
	return true;
    }
  
    bool VaneCalorimeter::isContainedSection(CLHEP::Hep3Vector const& start, CLHEP::Hep3Vector const& stop) const 
    {   
	for (unsigned int ivane=0;ivane<_nSections;++ivane) 
	   if (isInsideSection(ivane,start) && isInsideSection(ivane,stop)) return true;
	return false;
    }

   
    int VaneCalorimeter::crystalIdxFromPosition(const CLHEP::Hep3Vector &pos) const 
    {   
        int offset(0);
        for (unsigned int ivane=0;ivane<_nSections;++ivane)
	{
           if ( isInsideSection(ivane,pos) )
	   {
                 CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);
                 return offset + vane(ivane).idxFromPosition(posInSection.x(),posInSection.y());
  	   }
           offset +=vane(ivane).nCrystals();
        }       	  
        return -1;
    }

    std::vector<int> VaneCalorimeter::neighborsByLevel(int crystalId, int level) const 
    {

	int iv = _fullCrystalList.at(crystalId)->sectionId();

	int offset(0);
	for (int i=0;i<iv;++i) offset += vane(i).nCrystals();

	std::vector<int> list = vane(iv).findLocalNeighbors( _fullCrystalList.at(crystalId)->localId() ,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  

	return list;
    }

    void VaneCalorimeter::print(std::ostream &os) const 
    {
      os<<"Vane calorimeter "<<std::endl;
      os<<"Number of vanes :"<< _nSections<<std::endl;
      for (unsigned int ivane=0;ivane<_nSections;++ivane) vane(ivane).print(os);
        
    }



}
