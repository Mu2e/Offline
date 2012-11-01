// Add random noise hits with the specified probability.
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Digitization_PixelNoise_hh
#define ExtinctionMonitorFNAL_Digitization_PixelNoise_hh

#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

namespace mu2e {
  namespace ExtMonFNAL { class ExtMon; }

  class PixelNoise {
  public:

    PixelNoise(art::RandomNumberGenerator::base_engine_t& rng,
               const ExtMonFNAL::ExtMon **em,/*ugly workaroud  for Geometry not available at module ctr*/
               double noisePerPixelPerBC,
               int clockMin,
               int clockMax);

    // Adds noise hits to the collection.  Noise hits can be
    // identified by the lack of associated truth.  No need to have a
    // truth collection ptr arg.
    void add(ExtMonFNALRawHitCollection *hits);

  private:
    CLHEP::RandPoissonQ poisson_;
    CLHEP::RandFlat flat_;
    const ExtMonFNAL::ExtMon **extmon_;
    double noisePerPixelPerBC_;
    int clockMin_;
    int clockMax_;
  };

} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Digitization_PixelNoise_hh*/
