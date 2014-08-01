// $Id: DiskCrystalPosUtil.cc,v 1.2 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Utility to related the crystal indices in the disk to that of the Hexagonal map 
// it is needed because there is a hole in the disk, creating a mismatch between the numbering

// C++ includes
#include <iostream>
#include <map>

// Mu2e includes
#include "CalorimeterGeom/inc/DiskCrystalPosUtil.hh"



namespace mu2e {

      
      
      
      void DiskCrystalPosUtil::Fill(int iCrystal) {_mapToCrystal.push_back(iCrystal);}

      void DiskCrystalPosUtil::Fill(int iCrystal, int iMap, int l, int k) 
      {	    
	  _mapToCrystal.push_back(iMap);
	  _crystalToMap.push_back(iCrystal);
	  _crystalToL.push_back(l);
	  _crystalToK.push_back(k);
	  _LKToCrystal.insert( std::pair<int,int>(10000*l+k,iMap) );
      }

      int DiskCrystalPosUtil::mapToCrystal(int iMap)     const {return _mapToCrystal.at(iMap);}      
      int DiskCrystalPosUtil::crystalToMap(int iCrystal) const {return _crystalToMap.at(iCrystal);}

      int DiskCrystalPosUtil::crystalToL(int iCrystal)   const {return _crystalToL.at(iCrystal);}
      int DiskCrystalPosUtil::crystalToK(int iCrystal)   const {return _crystalToK.at(iCrystal);}      
      int DiskCrystalPosUtil::LKToCrystal(int l, int k)  const 
      { 
         std::map<int,int>::const_iterator ifind = _LKToCrystal.find(10000*l+k); 
	 return  ( ifind != _LKToCrystal.end() ) ? ifind->second : -1; 
      }

}


             
