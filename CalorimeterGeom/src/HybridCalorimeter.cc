//
// Geometry and identifier info about the Calorimeter.
//
// Original author B. Echenard
//

// C++ includes
#include <iostream>
#include <algorithm>

// Mu2e includes
#include "CalorimeterGeom/inc/HybridCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Barrel.hh"

//other includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {



   bool HybridCalorimeter::isInsideDisk(CLHEP::Hep3Vector const& pos) const 
   {   
	
	CLHEP::Hep3Vector posInSection = toSectionFrame(0, pos);
        
        double zlim    = _crystalHalfLength + _wrapperThickness + _caseThickness + _roHalfThickness + 0.5;         
	double rinlim  = _diskInnerRadius - _caseThickness - 0.5;
	double routlim = _diskOuterRadius + 0.5 + _caseThickness;
	
	if (posInSection.z()< -zlim || posInSection.z() > zlim ) return false;
	if (posInSection.perp() <  rinlim || posInSection.perp()  > routlim ) return false;
	return true;
    }


  bool HybridCalorimeter::isInsideBarrel(CLHEP::Hep3Vector const& pos) const 
   {   
	
	CLHEP::Hep3Vector posInSection = toSectionFrame(1, pos);
        
        double zlim    = (_crystalHalfTrans + _wrapperThickness + _caseThickness)*_nWheels;         
	double rinlim  = _barrelInnerRadius - _caseThickness - 0.5;
	double routlim =  _barrelInnerRadius + 2.0*(_crystalHalfLength + _wrapperThickness + _caseThickness + _roHalfThickness);
	
	if (posInSection.z()< -zlim || posInSection.z() > zlim ) return false;
	if (posInSection.perp() <  rinlim || posInSection.perp()  > routlim ) return false;
	return true;
    }

    bool HybridCalorimeter::isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const 
    {   
      if(isInsideDisk(pos) || isInsideBarrel(pos)) return true;
      return false;
    }
 
  bool HybridCalorimeter::isInsideSection(int iSec, CLHEP::Hep3Vector const& pos) const{

    if(iSec==0) return isInsideDisk(pos);
    if(iSec==1) return isInsideBarrel(pos);
    return false;
  }
 


    int HybridCalorimeter::crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const 
    {   
      if ( isInsideDisk(pos) ) {
	CLHEP::Hep3Vector posInSection = toSectionFrame(0, pos);
	return  disk().idxFromPosition(posInSection.x(),posInSection.y());
      } 
      
      if ( isInsideBarrel(pos) ) {
	CLHEP::Hep3Vector posInSection = toSectionFrame(1, pos);
	return  barrel().idxFromPosition(posInSection);
      } 
      
      return -1;
    }

    std::vector<int> HybridCalorimeter::neighbors(int crystalId, int level)  const
    {

	int iv = caloSectionId(crystalId);
	int ic = localCrystalId(crystalId);

	int offset(0);
	if(iv==1) offset += disk().nCrystals();

	std::vector<int> list;
	if(iv==0){
	  list = disk().neighbors(ic,level);
	}else if(iv==1){
	  list = barrel().neighbors(ic,level);
	}
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  

	return list;
    }



   double HybridCalorimeter::crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const
   {   
	CLHEP::Hep3Vector posInSection = toCrystalFrame(crystalId, pos);
	return posInSection.z();
   }





}
