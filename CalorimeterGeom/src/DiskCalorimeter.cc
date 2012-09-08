//
// Geometry and identifier info about the Calorimeter.
//
// Original author B. Echenard
//

// C++ includes
#include <iostream>
#include <cmath>
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/HexPositionMap.hh"
#include "CalorimeterGeom/inc/HexPosition.hh"
#include "CalorimeterGeom/inc/Disk.hh"




namespace mu2e {

   DiskCalorimeter::DiskCalorimeter()
   {
   }



   //get total number of readouts
   unsigned int DiskCalorimeter::nRO(void) const 
   {
      unsigned total(0);
      for (unsigned int i=0;i<_nDisks;++i) total += _disks[i].getCrystalMap().nCrystals();
      return total*_nROPerCrystal;
   }


   int DiskCalorimeter::getDiskId(int id) const
   {     
     for (unsigned int i=0;i<_nDisks;++i) {
       if (id < _disks[i].getCrystalMap().nCrystals()) return i;
       id -= _disks[i].getCrystalMap().nCrystals();
     }
     return _nDisks;
   }

   int DiskCalorimeter::getCrystalIdInMap(int id) const
   {     
     for (unsigned int i=0;i<_nDisks;++i) {
       if (id < _disks[i].getCrystalMap().nCrystals()) return id;
       id -= _disks[i].getCrystalMap().nCrystals();
     }
     return id;
   }



   // Get crystal origin (center) in Mu2e coordinates
   CLHEP::Hep3Vector DiskCalorimeter::getCrystalOriginById(int id) const 
   {
       int idisk = getDiskId(id);
       int ic    = getCrystalIdInMap(id);

              
       // Crystal center in disk coordinates
       CLHEP::Hep2Vector xypos = _disks[idisk].getCrystalMap().getCrystalPos(ic).XY();
       CLHEP::Hep3Vector plocal(xypos.x(),xypos.y(),-1.0*_readOutHalfThickness);
 
 
        //must apply the hack to get into Mu2e coordinates
       CLHEP::Hep3Vector  hack(-3904,0,0);
       CLHEP::Hep3Vector diskOrig = _disks[idisk].getOrigin()+hack;
      
       return diskOrig + (_disks[idisk].getRotation().inverse())*plocal; 
   }




   // Convert coordinates from Mu2e frame to local frame of crystal (taken as center of crystal), identified by id
   CLHEP::Hep3Vector DiskCalorimeter::toCrystalFrame(int roid, CLHEP::Hep3Vector const& pos) const 
   {   

       int id    = getCrystalByRO(roid);       
       int idisk = getDiskId(id);
       int ic    = getCrystalIdInMap(id);
       
       // Crystal center in disk coordinates
       CLHEP::Hep2Vector xypos = _disks[idisk].getCrystalMap().getCrystalPos(ic).XY();
       CLHEP::Hep3Vector plocal(xypos.x(),xypos.y(),-1.0*_readOutHalfThickness);
       
       //if you want to coordinates w.r.t the front face of the crystal, use
       //CLHEP::Hep3Vector plocal(xy.first,xy.second,-1.0*_readOutHalfThickness-_crystalDepth/2.0);
        
       
       //must apply the hack to get into Mu2e coordinates
       CLHEP::Hep3Vector  hack(-3904,0,0);
       CLHEP::Hep3Vector diskOrig = _disks[idisk].getOrigin()+hack;


       return (_disks[idisk].getRotation())*(pos-diskOrig)-plocal;  
   }




   CLHEP::Hep3Vector DiskCalorimeter::toDiskFrame(int idisk, CLHEP::Hep3Vector const& pos) const 
   {   
       CLHEP::Hep3Vector  hack(-3904,0,0);
       CLHEP::Hep3Vector diskOrig = _disks[idisk].getOrigin()+hack;
       return (_disks[idisk].getRotation())*(pos-diskOrig);
   }

   CLHEP::Hep3Vector DiskCalorimeter::fromDiskFrame(int idisk, CLHEP::Hep3Vector const& pos) const 
   {   
       CLHEP::Hep3Vector  hack(-3904,0,0);
       CLHEP::Hep3Vector diskOrig = _disks[idisk].getOrigin()+hack;
             
       return diskOrig + (_disks[idisk].getRotation().inverse())*pos;
   }




} // namespace mu2e
