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

double DiskCalorimeter::envelopeRmin() const {
	int ik=0;
	double radius_tmp[_nDisk];
	for (unsigned int i=0;i<_nDisk;++i) {
		radius_tmp[i]= _disks.at(i).innerRadius();
		if (i>0 && radius_tmp[i] < radius_tmp[i-1]) ik=i;
	}
	return _disks.at(ik).innerRadius() - _disks.at(ik).thickness() - 50.0 ;
}

double DiskCalorimeter::envelopeRmax() const {
	int ik=0;
	double radius_tmp[_nDisk];
	for (unsigned int i=0;i<_nDisk;++i) {
		radius_tmp[i]= _disks.at(i).innerRadius();
		if (i>0 && radius_tmp[i] > radius_tmp[i-1]) ik=i;
	}
	return _disks.at(ik).outerRadius() + _disks.at(ik).thickness() + 50.0;
}

double DiskCalorimeter::envelopeHalfLength() const {
	double diskTotThick = _crystalDepth + 2.0*_readOutHalfThickness + 2.0*_wrapperThickness+2.0*_diskThickness;
	return ( 2.0*diskTotThick + _diskSeparation.back() )/2.0 + 800.0;// 5.0;
}










    CLHEP::Hep3Vector DiskCalorimeter::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   

       const Disk& thisDisk = disk( caloSectionId(CrystalId) );
       int ic               = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisDisk.crystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       crysLocalPos += shift;

CLHEP::Hep3Vector hack(0,0,-600);
crysLocalPos += hack;
               
       return (thisDisk.rotation())*(pos-thisDisk.origin())-crysLocalPos;  
    }


    CLHEP::Hep3Vector DiskCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const Disk& thisDisk = disk(sectionId);
CLHEP::Hep3Vector hack(0,0,-600);
       return (thisDisk.rotation())*(pos-thisDisk.origin()) -hack;
    }


    CLHEP::Hep3Vector DiskCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const Disk& thisDisk = disk(sectionId);
CLHEP::Hep3Vector hack(0,0,-600);
        return thisDisk.inverseRotation()*pos + thisDisk.origin() +hack;
    }


    CLHEP::Hep3Vector DiskCalorimeter::crystalOrigin(int CrystalId) const 
    {          
       const Disk& thisDisk = disk( caloSectionId(CrystalId) );
       int ic               = localCrystalId(CrystalId);
             
       CLHEP::Hep3Vector crysLocalPos = thisDisk.crystal(ic).position(); 
       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       crysLocalPos += shift;

CLHEP::Hep3Vector hack(0,0,-600);
crysLocalPos -= hack;
                   
       return thisDisk.origin() + thisDisk.inverseRotation()*crysLocalPos; 
    }


    CLHEP::Hep3Vector DiskCalorimeter::localCrystalOrigin(int CrystalId) const 
    {          
       const Disk& thisDisk = disk( caloSectionId(CrystalId) );
       int ic               = localCrystalId(CrystalId);
             
       CLHEP::Hep3Vector crysLocalPos = thisDisk.crystal(ic).position(); 
       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       //CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       //crysLocalPos += shift;
                   
       return crysLocalPos; 
    }

    CLHEP::Hep3Vector DiskCalorimeter::crystalAxis(int CrystalId) const 
    {
       const Disk& thisDisk = disk( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector vlocal(0,0,1);
       return thisDisk.inverseRotation()*vlocal;
    }









    bool DiskCalorimeter::isInsideDisk(int idisk, CLHEP::Hep3Vector const& pos) const 
    {   
	
	CLHEP::Hep3Vector posInSection = toSectionFrame(idisk, pos);
        
        double zlim    = _crystalDepth/2.0 + _wrapperThickness + _diskThickness + _readOutHalfThickness + 0.5;         
	double rinlim  = _diskInnerRadius[idisk] - _diskThickness - 0.5;
	double routlim = _diskOuterRadius[idisk]+0.5 + _diskThickness;
	
	if (posInSection.z()< -zlim || posInSection.z() > zlim ) return false;
	if (posInSection.perp() <  rinlim || posInSection.perp()  > routlim ) return false;
	return true;
    }

    int DiskCalorimeter::crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const 
    {   
       int offset(0);
       for (unsigned int idisk=0;idisk<_nDisk;++idisk) {
 	  if ( isInsideDisk(idisk,pos) ) {
	   	CLHEP::Hep3Vector posInSection = toSectionFrame(idisk, pos);
		return offset + disk(idisk).idxFromPosition(posInSection.x(),posInSection.y());
	  } 
	  offset +=disk(idisk).nCrystals();
       }       	  
       return -1;
    }


    bool DiskCalorimeter::isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const 
    {   
       for (unsigned int idisk=0;idisk<_nDisk;++idisk) if (isInsideDisk(idisk, pos)) return true;
       return false;    
    }
 

 


    std::vector<int> DiskCalorimeter::neighbors(int crystalId, int level)  const
    {

       int iv = caloSectionId(crystalId);
       int ic = localCrystalId(crystalId);
       
       int offset(0);
       for (int i=0;i<iv;++i) offset +=disk(i).nCrystals();
              
       std::vector<int> list = disk(iv).neighbors(ic,level);
       transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  

       return list;
    }





    //get total number of readouts
    unsigned int DiskCalorimeter::nRO(void) const 
    {
       unsigned total(0);
       for (unsigned int i=0;i<_nDisk;++i) total += disk(i).nCrystals();
       return total*_nROPerCrystal;
    }

    unsigned int DiskCalorimeter::nCrystal(void) const 
    {
       unsigned total(0);
       for (unsigned int i=0;i<_nDisk;++i) total += disk(i).nCrystals();
       return total;
    }
    
    


    int DiskCalorimeter::caloSectionId(int crystalId) const
    {          
      for (unsigned int i=0;i<_nDisk;++i) 
      {
        if (crystalId < disk(i).nCrystals()) return i;
        crystalId -= disk(i).nCrystals();
      }
      return _nDisk;
    }
    
    
    int DiskCalorimeter::localCrystalId(int crystalId) const
    {     
      for (unsigned int i=0;i<_nDisk;++i) 
      {
        if (crystalId < disk(i).nCrystals()) return crystalId;
        crystalId -= disk(i).nCrystals();
      }
      return crystalId;
    }











}
