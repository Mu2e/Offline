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

double DiskCalorimeter::envelopeRmin(void) const {
	int ik=0;
	double radius_tmp[_nDisks];
	for (unsigned int i=0;i<_nDisks;++i) {
		radius_tmp[i]= _disks.at(i).innerRadius();
		if (i>0 && radius_tmp[i] < radius_tmp[i-1]) ik=i;
	}
	return _disks.at(ik).innerRadius() - _disks.at(ik).thickness() - 50.0 ;
}

double DiskCalorimeter::envelopeRmax(void) const {
	int ik=0;
	double radius_tmp[_nDisks];
	for (unsigned int i=0;i<_nDisks;++i) {
		radius_tmp[i]= _disks.at(i).innerRadius();
		if (i>0 && radius_tmp[i] > radius_tmp[i-1]) ik=i;
	}
	return _disks.at(ik).outerRadius() + _disks.at(ik).thickness() + 50.0;
}

double DiskCalorimeter::envelopeHalfLength(void) const {
	double diskTotThick = _crystalDepth + 2.0*_readOutHalfThickness + 2.0*_wrapperThickness+2.0*_diskThickness;
	return ( 2.0*diskTotThick + _diskSeparation.back() )/2.0 + 800.0;// 5.0;
}


    CLHEP::Hep3Vector DiskCalorimeter::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   

       const Disk& disk0 = disk( caloSectionId(CrystalId) );
       int ic           = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = disk0.crystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       crysLocalPos += shift;
               
       return (disk0.rotation())*(pos-disk0.origin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector DiskCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const Disk& disk0 = disk(sectionId);
       return (disk0.rotation())*(pos-disk0.origin());
    }

    CLHEP::Hep3Vector DiskCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const Disk& disk0 = disk(sectionId);
        return disk0.inverseRotation()*pos + disk0.origin();
    }

    CLHEP::Hep3Vector DiskCalorimeter::crystalOrigin(int CrystalId) const 
    {          
       const Disk& disk0 = disk( caloSectionId(CrystalId) );
       int ic           = localCrystalId(CrystalId);
             
       CLHEP::Hep3Vector crysLocalPos = disk0.crystal(ic).position(); 
       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       crysLocalPos += shift;
                  
       return disk0.origin() + disk0.inverseRotation()*crysLocalPos; 
    }

    CLHEP::Hep3Vector DiskCalorimeter::localCrystalOrigin(int CrystalId) const 
    {          
       const Disk& disk0 = disk( caloSectionId(CrystalId) );
       int ic           = localCrystalId(CrystalId);
             
       CLHEP::Hep3Vector crysLocalPos = disk0.crystal(ic).position(); 
       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       //CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       //crysLocalPos += shift;
                  
       return crysLocalPos; 
    }

    CLHEP::Hep3Vector DiskCalorimeter::crystalAxis(int CrystalId) const 
    {
       const Disk& disk0 = disk( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector vlocal(0,0,1);
       return disk0.inverseRotation()*vlocal;
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
       for (unsigned int i=0;i<_nDisks;++i) total += disk(i).nCrystals();
       return total*_nROPerCrystal;
    }

    unsigned int DiskCalorimeter::nCrystal(void) const 
    {
       unsigned total(0);
       for (unsigned int i=0;i<_nDisks;++i) total += disk(i).nCrystals();
       return total;
    }
    
    
    int DiskCalorimeter::caloSectionId(int crystalId) const
    {          
      for (unsigned int i=0;i<_nDisks;++i) {
        if (crystalId < disk(i).nCrystals()) return i;
        crystalId -= disk(i).nCrystals();
      }
      return _nDisks;
    }
    
    
    int DiskCalorimeter::localCrystalId(int crystalId) const
    {     
      for (unsigned int i=0;i<_nDisks;++i) {
        if (crystalId < disk(i).nCrystals()) return crystalId;
        crystalId -= disk(i).nCrystals();
      }
      return crystalId;
    }





}
