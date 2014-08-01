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
//C++ includes
#include <algorithm>

//mu2e includes
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"



namespace mu2e {





    bool VaneCalorimeter::isInsideVane(int ivane, CLHEP::Hep3Vector const& pos) const 
    {   
	CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);

	double xlim = _caloGeomInfo.crystalHalfLength() + _caloGeomInfo.wrapperThickness() + _caloGeomInfo.roHalfThickness()  + 0.5;
	double ylim = _nCrystalR*(_caloGeomInfo.crystalHalfTrans() +_caloGeomInfo.wrapperThickness()+_caloGeomInfo.shellThickness()) + 0.5;   
	double zlim = _nCrystalZ*(_caloGeomInfo.crystalHalfTrans() +_caloGeomInfo.wrapperThickness()+_caloGeomInfo.shellThickness()) + 0.5;   

	if (posInSection.x() < -xlim || posInSection.x() > xlim ) return false;      
	if (posInSection.y() < -ylim || posInSection.y() > ylim ) return false;      
	if (posInSection.z() < -zlim || posInSection.z() > zlim ) return false;      

	return true;
    }
  
    bool VaneCalorimeter::isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const 
    {   
        for (unsigned int ivane=0;ivane<_nSections;++ivane) if (isInsideVane(ivane,pos)) return true;
        return false;    
    }

   
    int VaneCalorimeter::crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const 
    {   
        int offset(0);
        for (unsigned int ivane=0;ivane<_nSections;++ivane) {
           if ( isInsideVane(ivane,pos) ) {
                 CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);
                 return offset + vane(ivane).idxFromPosition(posInSection.y(),posInSection.z());
  	   }
           offset +=vane(ivane).nCrystals();
        }       	  
        return -1;
    }

    std::vector<int> VaneCalorimeter::neighborsByLevel(int crystalId, int level) const 
    {

	int iv = caloSectionId(crystalId);

	int offset(0);
	for (int i=0;i<iv;++i) offset += vane(i).nCrystals();

	std::vector<int> list = vane(iv).findLocalNeighbors( _fullCrystalList.at(crystalId)->localId() ,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  

	return list;

    }


   double VaneCalorimeter::crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const
   {   
	CLHEP::Hep3Vector posInSection = toCrystalFrame(crystalId, pos);
	return posInSection.x();
   }


}
