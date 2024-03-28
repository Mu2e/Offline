#ifndef CosmicRayShieldGeom_CRSScintillatorBar_hh
#define CosmicRayShieldGeom_CRSScintillatorBar_hh
//
// Representation of one Scintillator Bar in CosmicRayShield.
//
//
// Original author KLG; somewhat based on Rob Kutschke's Straw
//

#include <vector>
#include <string>
#include <memory>

#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorBarId.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e
{

  // Forward declarations.
  class CosmicRayShield;

  class CRSScintillatorBar
  {

    friend class CRSScintillatorLayer;
    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShield;
//    friend class CosmicRayShieldMaker;

    CRSScintillatorBar();

    public:

    CRSScintillatorBar(CRSScintillatorBarIndex const &index,
                       CRSScintillatorBarId const &id,
                       CLHEP::Hep3Vector const &position,
                       const std::shared_ptr<CRSScintillatorBarDetail> detail);

    // Accept the compiler generated destructor, copy constructor and assignment operators
    CRSScintillatorBarIndex index() const {return _index;}
    CRSScintillatorBarId id() const {return _id;}

    // Formatted string embedding the id of the ScintillatorBar.
    std::string name( std::string const & base ) const;

    CLHEP::Hep3Vector const & getPosition() const {return _position;}
    std::vector<double> const & getHalfLengths() const {return _detail->getHalfLengths();}
    std::string const & getMaterialName() const {return _detail->getMaterialName();}

    double getHalfThickness() const { return _detail->getHalfThickness();}
    double getHalfWidth() const { return _detail->getHalfWidth();}
    double getHalfLength() const { return _detail->getHalfLength();}

    const CRSScintillatorBarDetail& getBarDetail() const
    {
      return *_detail;
    }

    CLHEP::Hep3Vector toWorld(const CLHEP::Hep3Vector &localPosition) const
    {
      return _detail->toWorld(localPosition,_position);
    }
    CLHEP::Hep3Vector toLocal(const CLHEP::Hep3Vector &worldPosition) const
    {
      return _detail->toLocal(worldPosition,_position);
    }
    CLHEP::Hep3Vector toLocalNormalized(const CLHEP::Hep3Vector &worldPosition) const
    {
      return _detail->toLocalNormalized(worldPosition,_position);
    }
    bool isInside(const CLHEP::Hep3Vector &worldPosition) const
    {
      return _detail->isInside(worldPosition,_position);
    }

    bool operator==(const CRSScintillatorBar &other) const
    {
      return _index == other.index();
    }
    bool operator>(const CRSScintillatorBar &other) const
    {
      return _index > other.index();
    }
    bool operator<(const CRSScintillatorBar &other) const
    {
      return _index < other.index();
    }

    private:

    // Index into the container of all ScintillatorBars.
    CRSScintillatorBarIndex _index;

    // Identifier within one layer, module, shield, and shield ID
    CRSScintillatorBarId _id;

    // Mid-point of the ScintillatorBar, in Mu2e
    CLHEP::Hep3Vector _position;

    // Detailed description of a bar
    std::shared_ptr<CRSScintillatorBarDetail> _detail;


    /********************/
    // counter motherboard section

    public:

    CLHEP::Hep3Vector getCMBPosition(int side) const  //"side" can be 0 or 1 (i.e. the CMBs on the negative and positive side of the counter)
    {
      return _detail->getCMBPosition(side, _position);
    }
    std::vector<double> getCMBHalfLengths() const
    {
      return _detail->getCMBHalfLengths();
    }
    std::string const & getCMBMaterialName() const {return _detail->getCMBMaterialName();}
    bool hasCMB(int side) const {return _detail->hasCMB(side);}

    /********************/
    // SiPM section

    public:

    CLHEP::Hep3Vector getSiPMPosition(int SiPMNumber) const
    {
      return _detail->getSiPMPosition(SiPMNumber, _position);
    }
    CLHEP::Hep3Vector getSiPMPosition(int fiberNumber, int side) const
    {
      return _detail->getSiPMPosition(fiberNumber, side, _position);
    }

  };

}  //namespace mu2e

#endif /* CosmicRayShieldGeom_CRSScintillatorBar_hh */
