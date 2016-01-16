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


    BaseCalorimeter::BaseCalorimeter() : 
      _caloType(Calorimeter::CaloType::none),
      _nSections(0),
      _sections(),
      _origin(),
      _trackerCenter(),
      _caloGeomInfo(),
      _fullCrystalList()  
    {}


    CLHEP::Hep3Vector BaseCalorimeter::toCrystalFrame(int crystalId, const CLHEP::Hep3Vector &pos) const 
    {   
        const CaloSection &thisSection = section( _fullCrystalList.at(crystalId)->sectionId());
        CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(crystalId)->localPosition();
        return thisSection.rotation()*(pos-thisSection.origin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector BaseCalorimeter::toSectionFrame(int sectionId, const CLHEP::Hep3Vector &pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return (thisSection.rotation())*(pos-thisSection.origin());
    }

    CLHEP::Hep3Vector BaseCalorimeter::toSectionFrameFF(int sectionId, const CLHEP::Hep3Vector &pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return (thisSection.rotation())*(pos-thisSection.origin()) - thisSection.originToCrystalOrigin();
    }

    CLHEP::Hep3Vector BaseCalorimeter::toTrackerFrame(const CLHEP::Hep3Vector &pos) const 
    {   
        return pos - _trackerCenter;
    }



    CLHEP::Hep3Vector BaseCalorimeter::fromCrystalFrame(int crystalId, const CLHEP::Hep3Vector &pos) const 
    {   
        const CaloSection& thisSection = section( _fullCrystalList.at(crystalId)->sectionId() );
        CLHEP::Hep3Vector crysLocalPos = _fullCrystalList.at(crystalId)->localPosition();
        return thisSection.inverseRotation()*(pos+crysLocalPos) + thisSection.origin();  
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromSectionFrame(int sectionId, const CLHEP::Hep3Vector &pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*pos + thisSection.origin();
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromSectionFrameFF(int sectionId, const CLHEP::Hep3Vector &pos) const 
    {   
        const CaloSection& thisSection = section(sectionId);
        return thisSection.inverseRotation()*(pos + thisSection.originToCrystalOrigin()) + thisSection.origin();
    }

    CLHEP::Hep3Vector BaseCalorimeter::fromTrackerFrame(const CLHEP::Hep3Vector &pos) const 
    {   
        return pos + _trackerCenter;
    }







}


