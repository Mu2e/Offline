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


    CLHEP::Hep3Vector BaseCalorimeter::toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section( _fullCrystalList.at(crystalId)->sectionId());
       CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(crystalId)->localPosition();
       return thisSection.rotation()*(pos-thisSection.origin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector BaseCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section(sectionId);
       return (thisSection.rotation())*(pos-thisSection.origin());
    }

    CLHEP::Hep3Vector BaseCalorimeter::toSectionFrameFF(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section(sectionId);
       return (thisSection.rotation())*(pos-thisSection.origin()) - thisSection.originToCrystalOrigin();
    }

    CLHEP::Hep3Vector BaseCalorimeter::toTrackerFrame(CLHEP::Hep3Vector const& pos) const 
    {   
        return pos - _trackerCenter;
    }



    CLHEP::Hep3Vector BaseCalorimeter::fromCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const 
    {   
       const CaloSection& thisSection = section( _fullCrystalList.at(crystalId)->sectionId() );
       CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(crystalId)->localPosition();
       return thisSection.inverseRotation()*(pos+crysLocalPos) + thisSection.origin();  
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*pos + thisSection.origin();
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromSectionFrameFF(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*(pos + thisSection.originToCrystalOrigin()) + thisSection.origin();
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromTrackerFrame(CLHEP::Hep3Vector const& pos) const 
    {   
        return pos + _trackerCenter;
    }




    CLHEP::Hep3Vector BaseCalorimeter::crystalOrigin(int crystalId) const 
    {          
       // -- this is the original code. For efficiency purpose, I precalculate the position in the mu2e frame 
       // -- and return it instead of recalculating it all the time
       //const CaloSection& thisSection = section( _fullCrystalList.at(crystalId)->sectionId());
       //return thisSection.origin() + thisSection.inverseRotation()*_fullCrystalList.at(crystalId)->localPosition();        
       return _fullCrystalList.at(crystalId)->position();    
    }

    CLHEP::Hep3Vector BaseCalorimeter::crystalOriginInSection(int crystalId) const 
    {          
       return _fullCrystalList.at(crystalId)->localPosition(); 
    }


    CLHEP::Hep3Vector BaseCalorimeter::crystalAxis(int crystalId) const 
    {
       CLHEP::Hep3Vector vlocal(0,0,1);
       const CaloSection& thisSection = section(_fullCrystalList.at(crystalId)->sectionId());
       return thisSection.inverseRotation()*vlocal;
    }





}


