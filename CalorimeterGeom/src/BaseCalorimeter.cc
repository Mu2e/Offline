//
// Geometry and identifier info about the Calorimeter.
//
// Original author B. Echenard
//

// C++ includes
#include <iostream>
#include <algorithm>

// Mu2e includes
#include "CalorimeterGeom/inc/BaseCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "cetlib/exception.h"

//other includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {


    int BaseCalorimeter::caloSectionId(int crystalId) const
    {          
        if ( crystalId >= _nCrystalTot) throw cet::exception("Calorimeter") << "Tried to access crystal Id "<<crystalId<<" but maximum is "<<_nCrystalTot-1<<"! \n";      

        for (unsigned int i=0;i<_nSections;++i) 
        {
           if (crystalId < section(i).nCrystals()) return i;
           crystalId -= section(i).nCrystals();
        }
        return _nSections;
    }

    int BaseCalorimeter::localCrystalId(int crystalId) const
    {     
        if ( crystalId >= _nCrystalTot) throw cet::exception("Calorimeter") << "Tried to access crystal Id "<<crystalId<<" but maximum is "<<_nCrystalTot-1<<"! \n";      
        
	for (unsigned int i=0;i<_nSections;++i) 
        {
          if (crystalId < section(i).nCrystals()) return crystalId;
          crystalId -= section(i).nCrystals();
        }
        return crystalId;
    }





    unsigned int BaseCalorimeter::nRO(void) const 
    {
        unsigned total(0);
        for (unsigned int i=0;i<_nSections;++i) total += section(i).nCrystals();
        return total*_nROPerCrystal;
    }

    unsigned int BaseCalorimeter::nCrystal(void) const 
    {
        unsigned total(0);
        for (unsigned int i=0;i<_nSections;++i) total += section(i).nCrystals();
        return total;
    }







    CLHEP::Hep3Vector BaseCalorimeter::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   

       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       int ic                         = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisSection.crystal(ic).position();
       crysLocalPos += _crystalShift;

       return (thisSection.rotation())*(pos-thisSection.origin())-crysLocalPos;  
    }


    CLHEP::Hep3Vector BaseCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section(sectionId);
       return (thisSection.rotation())*(pos-thisSection.origin());
    }


    CLHEP::Hep3Vector BaseCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*pos + thisSection.origin();
    }


    CLHEP::Hep3Vector BaseCalorimeter::crystalOrigin(int CrystalId) const 
    {          
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       int ic                         = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisSection.crystal(ic).position();
       crysLocalPos += _crystalShift;

       return thisSection.origin() + thisSection.inverseRotation()*crysLocalPos; 
    }


    CLHEP::Hep3Vector BaseCalorimeter::localCrystalOrigin(int CrystalId) const 
    {          
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       int ic                         = localCrystalId(CrystalId);
       
       CLHEP::Hep3Vector crysLocalPos = thisSection.crystal(ic).position();
       crysLocalPos += _crystalShift;
                   
       return crysLocalPos; 
    }


    CLHEP::Hep3Vector BaseCalorimeter::crystalAxis(int CrystalId) const 
    {
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector vlocal(0,0,1);
       return thisSection.inverseRotation()*vlocal;
    }





}
