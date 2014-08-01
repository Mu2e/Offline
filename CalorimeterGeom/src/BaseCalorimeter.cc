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



    CLHEP::Hep3Vector BaseCalorimeter::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(CrystalId)->localPosition()+ thisSection.crystalShift();

       return thisSection.rotation()*(pos-thisSection.origin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector BaseCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section(sectionId);
       return (thisSection.rotation())*(pos-thisSection.origin());
    }


    CLHEP::Hep3Vector BaseCalorimeter::fromCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(CrystalId)->localPosition()+ thisSection.crystalShift();

       return thisSection.inverseRotation()*(pos+crysLocalPos) + thisSection.origin();  
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*pos + thisSection.origin();
    }


    CLHEP::Hep3Vector BaseCalorimeter::crystalOrigin(int CrystalId) const 
    {          
       // -- this is the original code. For efficiency purpose, I precalculate the position in the mu2e frame 
       // -- and return it instead of recalculating it all the time
       //const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       //CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(CrystalId)->localPosition() + thisSection.crystalShift();
       //return thisSection.origin() + thisSection.inverseRotation()*crysLocalPos; 
       
       return _fullCrystalList.at(CrystalId)->position();    
    }


    CLHEP::Hep3Vector BaseCalorimeter::crystalOriginInSection(int CrystalId) const 
    {          
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(CrystalId)->localPosition()+ thisSection.crystalShift();
                   
       return crysLocalPos; 
    }

    CLHEP::Hep3Vector BaseCalorimeter::crystalAxis(int CrystalId) const 
    {
       const CaloSection& thisSection = section( caloSectionId(CrystalId) );
       CLHEP::Hep3Vector vlocal(0,0,1);
       return thisSection.inverseRotation()*vlocal;
    }





}


