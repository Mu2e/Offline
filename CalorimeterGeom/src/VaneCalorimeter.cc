//
// Geometry and identifier info about the VaneCalorimeter.
//
//
// $Id: VaneCalorimeter.cc,v 1.6 2013/05/28 22:11:24 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/28 22:11:24 $
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

	double xlim = _crystalHL + _wrapperThickness + _roHalfThickness + 0.5;
	double ylim = _nCrystalR*(_crystalHW+_wrapperThickness+_shellThickness) + 0.5;   
	double zlim = _nCrystalZ*(_crystalHW+_wrapperThickness+_shellThickness) + 0.5;   

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

    std::vector<int> VaneCalorimeter::neighbors(int CrystalId, int level) const 
    {

	int iv = caloSectionId(CrystalId);
	int ic = localCrystalId(CrystalId);

        int offset(0);
        for (int i=0;i<iv;++i) offset += vane(i).nCrystals();

	std::vector<int> list = vane(iv).neighbors(ic,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  
	return list;
    }


   double VaneCalorimeter::crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const
   {   
	CLHEP::Hep3Vector posInSection = toCrystalFrame(crystalId, pos);
	return posInSection.x();
   }


}
