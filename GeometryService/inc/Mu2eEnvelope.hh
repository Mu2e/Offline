// This define a box, in the Mu2e coordinate system, which contains
// all the pieces of "real" mu2e geometry, such as detectors,
// building, and shielding dirt.
//
// Andrei Gaponenko, 2012

#ifndef MU2EENVELOPE_HH
#define MU2EENVELOPE_HH

#include <ostream>

#include "Mu2eInterfaces/inc/Detector.hh"
#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class GeometryService;
  class Mu2eHall;
  class SimpleConfig;
  class ProtonBeamDump;
  class ExtMonFNALBuilding;

  class Mu2eEnvelope : virtual public Detector {
  public:

    double xmin() const { return xmin_; }
    double xmax() const { return xmax_; }
    double ymin() const { return ymin_; }
    double ymax() const { return ymax_; }
    double zmin() const { return zmin_; }
    double zmax() const { return zmax_; }

    //----------------------------------------------------------------
  private:
    // Private ctr: the class should be only obtained via GeometryService
    friend class GeometryService;
    Mu2eEnvelope(const Mu2eHall&, const SimpleConfig&);

    template<class T> friend class art::Wrapper; // Needed for persistency
    Mu2eEnvelope(); // Needed for persistency

    double xmin_;
    double xmax_;
    double ymin_;
    double ymax_;
    double zmin_;
    double zmax_;

  };

  std::ostream& operator<<(std::ostream& os, const Mu2eEnvelope& env);

}

#endif/*MU2EENVELOPE_HH*/
