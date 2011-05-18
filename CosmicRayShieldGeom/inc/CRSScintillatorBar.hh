#ifndef CosmicRayShieldGeom_CRSScintillatorBar_hh
#define CosmicRayShieldGeom_CRSScintillatorBar_hh
//
// Representation of one Scintillator Bar in CosmicRayShield.
//
// $Id: CRSScintillatorBar.hh,v 1.4 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
//
// Original author KLG; somewhat based on Rob Kutschke's Straw
//

#include <vector>
#include <string>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarId.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarIndex.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"

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
                       CLHEP::Hep3Vector       const & localOffset,
                       std::vector<double>     const & globalRotationAngles,
                       CLHEP::Hep3Vector       const & globalOffset
                       );

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const CRSScintillatorBarId& Id() const { return _id;}
    CRSScintillatorBarIndex Index() const { return _index;}

    // Formatted string embedding the id of the ScintillatorBar.
    std::string name( std::string const & base ) const;

    CLHEP::Hep3Vector const & getLocalOffset() const {return _localOffset;}

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
      if (_index == other.Index()) {
	return true;
      }
      else{
	return false;
      }
    }
    bool operator>(const CRSScintillatorBar other) const {
      if (_index > other.Index()) {
	return true;
      }
      else{
	return false;
      }
    }
   bool operator<(const CRSScintillatorBar other) const {
      if (_index < other.Index()) {
	return true;
      }
      else{
	return false;
      }
   }

  private:

    // Identifier
    CRSScintillatorBarId _id;

    // Index into the container of all ScintillatorBars.
    CRSScintillatorBarIndex _index;

    // Mid-point of the ScintillatorBar, a.k.a. localOffset
    CLHEP::Hep3Vector _localOffset;

    std::vector<double> _globalRotationAngles;
    // we may do rotation instead

    // Mid-point of the ScintillatorBar, in Mu2e
    CLHEP::Hep3Vector _globalOffset;

    // Detailed description of a bar
    mutable const CRSScintillatorBarDetail* _detail;
    // not needed for one detail int32_t _detailIndex;

    //     Nearest neighbours.// not filled out yet
    //     std::vector<CRSScintillatorBarId>    _nearestById;
    //     std::vector<CRSScintillatorBarIndex> _nearestByIndex;

  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorBar_hh */
