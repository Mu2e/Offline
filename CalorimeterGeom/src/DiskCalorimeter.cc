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



    CLHEP::Hep3Vector DiskCalorimeter::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   

       const Disk& disk = getDisk( getCaloSectionId(CrystalId) );
       int ic           = getLocalCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = disk.getCrystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       crysLocalPos += shift;
               
       return (disk.getRotation())*(pos-disk.getOrigin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector DiskCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const Disk& disk = getDisk(sectionId);
       return (disk.getRotation())*(pos-disk.getOrigin());
    }

    CLHEP::Hep3Vector DiskCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const Disk& disk = getDisk(sectionId);
        return disk.getInverseRotation()*pos + disk.getOrigin();
    }

    CLHEP::Hep3Vector DiskCalorimeter::getCrystalOrigin(int CrystalId) const 
    {          
       const Disk& disk = getDisk( getCaloSectionId(CrystalId) );
       int ic           = getLocalCrystalId(CrystalId);
             
       CLHEP::Hep3Vector crysLocalPos = disk.getCrystal(ic).position(); 
       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       crysLocalPos += shift;
                  
       return disk.getOrigin() + disk.getInverseRotation()*crysLocalPos; 
    }

    CLHEP::Hep3Vector DiskCalorimeter::getLocalCrystalOrigin(int CrystalId) const 
    {          
       const Disk& disk = getDisk( getCaloSectionId(CrystalId) );
       int ic           = getLocalCrystalId(CrystalId);
             
       CLHEP::Hep3Vector crysLocalPos = disk.getCrystal(ic).position(); 
       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       //CLHEP::Hep3Vector shift(0,0,-_crystalDepth/2.0);
       //crysLocalPos += shift;
                  
       return crysLocalPos; 
    }

    CLHEP::Hep3Vector DiskCalorimeter::getCrystalAxis(int CrystalId) const 
    {
       const Disk& disk = getDisk( getCaloSectionId(CrystalId) );
       CLHEP::Hep3Vector vlocal(0,0,1);
       return disk.getInverseRotation()*vlocal;
    }





    std::vector<int> DiskCalorimeter::getNeighbors(int crystalId, int level)  const
    {

       int iv = getCaloSectionId(crystalId);
       int ic = getLocalCrystalId(crystalId);
       
       int offset(0);
       for (int i=0;i<iv;++i) offset +=getDisk(i).nCrystals();
              
       std::vector<int> list = getDisk(iv).getNeighbors(ic,level);
       transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  

       return list;
    }



    //get total number of readouts
    unsigned int DiskCalorimeter::nRO(void) const 
    {
       unsigned total(0);
       for (unsigned int i=0;i<_nDisks;++i) total += getDisk(i).nCrystals();
       return total*_nROPerCrystal;
    }

    unsigned int DiskCalorimeter::nCrystal(void) const 
    {
       unsigned total(0);
       for (unsigned int i=0;i<_nDisks;++i) total += getDisk(i).nCrystals();
       return total;
    }
    
    
    int DiskCalorimeter::getCaloSectionId(int crystalId) const
    {          
      for (unsigned int i=0;i<_nDisks;++i) {
        if (crystalId < getDisk(i).nCrystals()) return i;
        crystalId -= getDisk(i).nCrystals();
      }
      return _nDisks;
    }
    
    
    int DiskCalorimeter::getLocalCrystalId(int crystalId) const
    {     
      for (unsigned int i=0;i<_nDisks;++i) {
        if (crystalId < getDisk(i).nCrystals()) return crystalId;
        crystalId -= getDisk(i).nCrystals();
      }
      return crystalId;
    }





}
