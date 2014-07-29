#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCendBoxes_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCendBoxes_hh

//
//
// Original author David Norvil Brown
//

// This class describes a set of "boxes" comprising the endcap neutron
// shield for the DS.  Created initially from the G4bl description of 
// of Geometry 13.

#include <vector>
#include <ostream>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class ExtNeutShieldCendBoxesMaker;

  class ExtNeutShieldCendBoxes : virtual public Detector {

  public:

    // Use a vector of Hep2Vectors for the corners of the shape to be extruded
    const std::vector<std::vector<double>>& getDimensions() const { return _dimensions; }
    const std::vector<std::string>& materialNames() const { return _materialNames; }
    const std::vector<CLHEP::Hep3Vector>& centersOfBoxes() const { return _centerPositions; }
    const bool hasHole( int i ) const { return _hasHole[i]; }
    const int holeIndex( int i ) const { return _holeIndexes[i]; }
    const double holeRadius( int i ) const { return _holeRadius[i]; }
    const CLHEP::Hep3Vector& holeLocation( int i ) const { return _holeLocations[i]; }
    const double holeHalfLength( int i ) const { return _holeHalfLength[i]; }

  private:

    friend class ExtNeutShieldCendBoxesMaker;

    // Private ctr: the class should only be constructed via ExtNeutShieldCendBoxes::ExtNeutShieldCendBoxesMaker.
//     ExtNeutShieldCendBoxes(const std::vector<std::vector<double>>& dims, 
//                            const std::vector<std::string>& mats, 
//                            const std::vector<CLHEP::Hep3Vector>& sites, 
//                            const std::vector<bool>& holey, 
//                            const std::vector<int>& whatHole, 
//                            const std::vector<CLHEP::Hep3Vector>& holeLocs,
//                            const std::vector<double>& rads, 
//                            const std::vector<double>& hlengs )
//       : _dimensions(dims),
//         _materialNames(mats),
// 	_centerPositions(sites),
// 	_hasHole(holey),
// 	_holeIndexes(whatHole),
//         _holeLocations(holeLocs),
// 	_holeRadius(rads),
// 	_holeHalfLength(hlengs)
//     { }

    // Or read back from persistent storage
    ExtNeutShieldCendBoxes();
    template<class T> friend class art::Wrapper;

    // Current description based on Geometry 13 from G4Beamline, adapted by
    // D. Norvil Brown, 


    std::vector<std::vector<double>> _dimensions;
    std::vector<std::string> _materialNames;
    std::vector<CLHEP::Hep3Vector> _centerPositions;
    std::vector<bool> _hasHole;
    std::vector<int> _holeIndexes;
    std::vector<CLHEP::Hep3Vector> _holeLocations;
    std::vector<double> _holeRadius;
    std::vector<double> _holeHalfLength;
  };

  std::ostream& operator<<(std::ostream& os, const ExtNeutShieldCendBoxes& enscb);

}

#endif/*ExternalNeutronShieldingGeom_ExtNeutShieldCendBoxes_hh*/
