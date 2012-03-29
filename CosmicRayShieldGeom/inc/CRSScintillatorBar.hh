#ifndef CosmicRayShieldGeom_CRSScintillatorBar_hh
#define CosmicRayShieldGeom_CRSScintillatorBar_hh
//
// Representation of one Scintillator Bar in CosmicRayShield.
//
// $Id: CRSScintillatorBar.hh,v 1.9 2012/03/29 22:59:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/03/29 22:59:13 $
//
// Original author KLG; somewhat based on Rob Kutschke's Straw
//

#include <vector>
#include <string>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarId.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward declarations.
  class CosmicRayShield;

  class CRSScintillatorBar{

    friend class CRSScintillatorLayer;
    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShield;
    friend class CosmicRayShieldMaker;

  public:

    CRSScintillatorBar();

    CRSScintillatorBar(
                       CRSScintillatorBarId    const & id,
                       CRSScintillatorBarIndex const & index
                       );

    CRSScintillatorBar(
                       CRSScintillatorBarId    const & id,
                       CRSScintillatorBarIndex const & index,
                       std::vector<double>     const & globalRotationAngles,
                       CLHEP::Hep3Vector       const & globalOffset
                       );

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const CRSScintillatorBarId& id() const { return _id;}
    CRSScintillatorBarIndex index() const { return _index;}

    // Formatted string embedding the id of the ScintillatorBar.
    std::string name( std::string const & base ) const;

    std::vector<double> const & getGlobalRotationAngles() const { return _globalRotationAngles;}

    CLHEP::Hep3Vector const & getGlobalOffset() const {return _globalOffset;}

    // CRSScintillatorBar Half lengths
    std::vector<double> const & getHalfLengths() {
      return _detail->getHalfLengths();
    }

    // CRSScintillatorBar Material Names
    std::vector<std::string> const & getMaterialNames() {
      return _detail->getMaterialNames();
    }

    // On readback from persistency, recursively recompute mutable members.
    //    void fillPointers ( const CosmicRayShield& cosmicRayShield ) const;

    bool operator==(const CRSScintillatorBar other) const {
      return _index == other.index();
    }
    bool operator>(const CRSScintillatorBar other) const {
      return _index > other.index();
    }
    bool operator<(const CRSScintillatorBar other) const {
      return _index < other.index();
   }

  private:

    // Identifier
    CRSScintillatorBarId _id;

    // Index into the container of all ScintillatorBars.
    CRSScintillatorBarIndex _index;

    std::vector<double> _globalRotationAngles;
    // we may do rotation instead

    // Mid-point of the ScintillatorBar, in Mu2e
    CLHEP::Hep3Vector _globalOffset;

    // Detailed description of a bar
    mutable const CRSScintillatorBarDetail* _detail;

    //     Nearest neighbours.// not filled out yet
    //     std::vector<CRSScintillatorBarId>    _nearestById;
    //     std::vector<CRSScintillatorBarIndex> _nearestByIndex;

  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorBar_hh */
