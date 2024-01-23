// Add random noise hits with the specified probability.
//
// Andrei Gaponenko, 2012


#include "Offline/ExtinctionMonitorFNAL/Digitization/inc/PixelNoise.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <map>

#include "cetlib_except/exception.h"

#include "Offline/DataProducts/inc/ExtMonFNALPixelId.hh"
#include "Offline/DataProducts/inc/ExtMonFNALPixelDenseId.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelIdConverter.hh"

#include "Offline/RecoDataProducts/inc/ExtMonFNALRawHit.hh"

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {


    namespace {

      typedef std::multimap<ExtMonFNALPixelDenseId,ExtMonFNALRawHit> NoiseMap;

      bool isHigh(const ExtMonFNALRawHit& hit, int clock) {
        return (hit.clock() <= clock) && (clock <= hit.clock() + int(hit.tot()));
      }

      void removeConflictingNoise(NoiseMap *nm,
                                  const ExtMonFNALRawHit& signalHit,
                                  const ExtMonFNALPixelIdConverter& conv)
      {
        std::vector<NoiseMap::iterator> toBeDeleted;

        typedef std::pair<NoiseMap::iterator, NoiseMap::iterator> Range;
        Range r = nm->equal_range(conv.densePixelNumber(signalHit.pixelId()));

        for(NoiseMap::iterator i = r.first; i != r.second; ++i) {
          if(isHigh(signalHit, i->second.clock())) {
            toBeDeleted.push_back(i);
            AGDEBUG("Found conflicting noise hit: "<<i->second<<" for signal hit "<<signalHit);
          }
        }

        for(unsigned i=0; i<toBeDeleted.size(); ++i) {
          nm->erase(toBeDeleted[i]);
        }

      }

    } // namespace{}


    //================================================================
    PixelNoise::PixelNoise(art::RandomNumberGenerator::base_engine_t& rng,
                           const ExtMon **em,
                           const int numBCs,
                           double noisePerPixelPerBC)
      : poisson_(rng)
      , flat_(rng)
      , extmon_(em)
      , numBCs_(numBCs)
      , noisePerPixelPerBC_(noisePerPixelPerBC)
    {}

    //================================================================
    void PixelNoise::add(ExtMonFNALRawHitCollection *hits) {
      ExtMonFNALPixelIdConverter conv(**extmon_);

      const double meanNoise = conv.totalNumberOfPixels() * numBCs_ * noisePerPixelPerBC_;

      int numNoise = poisson_.fire(meanNoise);

      AGDEBUG("Adding "<<numNoise<<" noise hits, mean noise = "<<meanNoise<<std::endl);

      // Generate the noise hits
      NoiseMap noise;
      for(int i=0; i<numNoise; ++i) {
        ExtMonFNALPixelDenseId pix(flat_.fireInt(conv.totalNumberOfPixels()));
        unsigned int clock = flat_.fireInt(numBCs_);
        unsigned int tot = 0; // noise hits are short
        noise.insert(std::make_pair(pix, ExtMonFNALRawHit(conv.pixelId(pix), clock, tot)));
      }

      // Retain only those noise hits that do not overlap with an existing hit
      for(ExtMonFNALRawHitCollection::const_iterator i=hits->begin(); i!=hits->end(); ++i) {
        removeConflictingNoise(&noise, *i, conv);
      }

      // Add surviving noise hits to the output
      for(NoiseMap::const_iterator i=noise.begin(); i!=noise.end(); ++i) {
        hits->push_back(i->second);
      }

      AGDEBUG("Final hits size = "<<hits->size());
    }

  } // namespace ExtMonFNAL
} // namespace mu2e
