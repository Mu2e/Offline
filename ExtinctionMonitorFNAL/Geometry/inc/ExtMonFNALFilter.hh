// A class to provide information to construct the collimator
// filter for the ExtMonFNAL detector:  entrance collimator,
// filter magnet, and exit collimator.
//
// Andrei Gaponenko, 2024

#ifndef EXTMONFNALFILTER_HH
#define EXTMONFNALFILTER_HH

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALCollimator.hh"

namespace mu2e {

  class ExtMonFNALBuilding;
  class ExtMonFNALFilterMaker;

  class ExtMonFNALFilter {
  public:

    double nominalMomentum() const { return nominalMomentum_; }

    const ExtMonFNALMagnet& magnet() const { return magnet_; }
    const ExtMonFNALCollimator& collimator1() const { return collimator1_; }
    const ExtMonFNALCollimator& collimator2() const { return collimator2_; }

    // convenience accessors
    CLHEP::Hep3Vector entranceInMu2e() const { return collimator1_.entranceInMu2e(); }
    CLHEP::Hep3Vector exitInMu2e() const { return collimator2_.exitInMu2e(); }

  private:
    friend class ExtMonFNALFilterMaker;
    // Private ctr: the class should be only obtained via the maker
    ExtMonFNALFilter() : nominalMomentum_{0.} {}
    // But the building may start with an uninitialized instance
    friend class ExtMonFNALBuilding;

    ExtMonFNALMagnet magnet_;
    ExtMonFNALCollimator collimator1_;
    ExtMonFNALCollimator collimator2_;
    double nominalMomentum_;
  };


} // namespace mu2e

#endif/*EXTMONFNALFILTER_HH*/
