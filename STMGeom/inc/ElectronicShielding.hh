#ifndef STMGeom_ElectronicShielding_hh
#define STMGeom_ElectronicShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

// Update: Yongyi Wu Oct 2026 -- a second constructor for the hand-stacked
// geometry (STM_v11): lab standard concrete blocks stacked on the hall
// floor and two silicon tile grids, every position resolved in Mu2e
// coordinates by STMMaker.

#include <string>
#include <vector>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class ElectronicShielding {
  public:
    // The earlier description, placed in constructSTM.cc relative to the
    // front shielding.
    ElectronicShielding(bool build,
    double SiGridX, double SiGridY, double SiGridZ,
    double SiXcenter, double SiYcenter, double SiZcenter,
    double ConcreteT, double GapToSi
    ):
      _build(build),
     _SiGridX(SiGridX), _SiGridY(SiGridY), _SiGridZ(SiGridZ),
     _SiXcenter(SiXcenter), _SiYcenter(SiYcenter), _SiZcenter(SiZcenter),
     _ConcreteT(ConcreteT), _GapToSi(GapToSi)
    {}

    // The hand-stacked description. All concrete blocks share one size
    // and all silicon tiles another, each given by its half lengths;
    // every block and every grid is given by its center in Mu2e
    // coordinates.
    ElectronicShielding(bool build,
    CLHEP::Hep3Vector const & concreteHalfLengths,
    std::vector<CLHEP::Hep3Vector> const & concreteCentersInMu2e,
    CLHEP::Hep3Vector const & siTileHalfLengths,
    int siNX, int siNY,
    std::vector<CLHEP::Hep3Vector> const & siGridCentersInMu2e
    ):
      _build(build),
     _concreteHalfLengths(concreteHalfLengths),
     _concreteCentersInMu2e(concreteCentersInMu2e),
     _siTileHalfLengths(siTileHalfLengths),
     _siNX(siNX), _siNY(siNY),
     _siGridCentersInMu2e(siGridCentersInMu2e)
    {}

    bool    build()                              const {return _build;}
    double  SiGridX()                            const {return _SiGridX;}
    double  SiGridY()                            const {return _SiGridY;}
    double  SiGridZ()                            const {return _SiGridZ;}
    double  SiXcenter()                          const {return _SiXcenter;}
    double  SiYcenter()                          const {return _SiYcenter;}
    double  SiZcenter()                          const {return _SiZcenter;}
    double  ConcreteT()                          const {return _ConcreteT;}
    double  GapToSi()                            const {return _GapToSi;}

    // ---- hand-stacked only; zero or empty for the earlier one ---------
    CLHEP::Hep3Vector const & concreteHalfLengths()   const {return _concreteHalfLengths;}
    std::vector<CLHEP::Hep3Vector> const & concreteCentersInMu2e() const {return _concreteCentersInMu2e;}
    CLHEP::Hep3Vector const & siTileHalfLengths()     const {return _siTileHalfLengths;}
    int     siNX()                               const {return _siNX;}
    int     siNY()                               const {return _siNY;}
    std::vector<CLHEP::Hep3Vector> const & siGridCentersInMu2e()   const {return _siGridCentersInMu2e;}


    ElectronicShielding() {}

  private:

    bool               _build;
    double             _SiGridX = 0.;
    double             _SiGridY = 0.;
    double             _SiGridZ = 0.;
    double             _SiXcenter = 0.;
    double             _SiYcenter = 0.;
    double             _SiZcenter = 0.;
    double             _ConcreteT = 0.;
    double             _GapToSi = 0.;

    CLHEP::Hep3Vector  _concreteHalfLengths;
    std::vector<CLHEP::Hep3Vector> _concreteCentersInMu2e;
    CLHEP::Hep3Vector  _siTileHalfLengths;
    int                _siNX = 0;
    int                _siNY = 0;
    std::vector<CLHEP::Hep3Vector> _siGridCentersInMu2e;

  };

}

#endif/*STMGeom_ElectronicShielding_hh*/
