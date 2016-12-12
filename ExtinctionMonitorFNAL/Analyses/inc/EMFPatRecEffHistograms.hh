// Track finding efficiency histograms.
//
// Andrei Gaponenko
//

#ifndef ExtinctionMonitorFNAL_Analyses_EMFPatRecEffHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFPatRecEffHistograms_hh

#include <string>

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "boost/noncopyable.hpp"

#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"

class TH1D;
class TH2D;
class TEfficiency;

namespace art { class Event; }
namespace art { class TFileDirectory; }

namespace mu2e {

  class SimParticle;

  namespace ExtMonFNAL {

    class ExtMon;

    class EMFPatRecEffHistograms :  private boost::noncopyable {
    public:

      explicit EMFPatRecEffHistograms(unsigned cutMinCommonClusters);

      // Book histograms in the subdirectory, given by the relativePath; that path is
      // relative to the root TFileDirectory for the current module.
      // The default is to use the current module's TFileDirectory
      void book(const ExtMon& extmon, const std::string& relativePath="");

      // Book histograms in the specified TFileDirectory.
      void book(const ExtMon& extmon, art::TFileDirectory& tfdir);

      // A helper class to use per-event FinMany object
      // instead of creating one per fill call.
      class Fillable {
        EMFPatRecEffHistograms *parent_;
        art::FindMany<ExtMonFNALTrkFit,ExtMonFNALTrkMatchInfo> trackFinder_;
        unsigned multiplicity_; // per-event variable
      public:

        Fillable(EMFPatRecEffHistograms *parent,
                 const std::vector<art::Ptr<SimParticle> >& particles,
                 const art::Event& event,
                 const art::InputTag& trkTruthTag,
                 unsigned multiplicity)
          : parent_(parent)
          , trackFinder_(particles, event, trkTruthTag)
          , multiplicity_(multiplicity)
        {}

        // Arg is an index in the particles collection given to the ctr.
        // Returns true for "efficient", false for inefficiencies
        bool fill(unsigned denominatorParticleIndex);
      };

      Fillable fillable(const std::vector<art::Ptr<SimParticle> >& particles,
                        const art::Event& event,
                        const art::InputTag& trkTruthTag,
                        unsigned multiplicity)  {
        return Fillable(this, particles, event, trkTruthTag, multiplicity);
      }

    private:
      friend class Fillable;
      unsigned cutMinCommonClusters_; // inclusive

      ExtMon *extmon_; // non-owning

      TH1D *hcommon_;
      TEfficiency* effMultiplicity_;
    };

  } // end namespace ExtMonFNAL
} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFPatRecEffHistograms_hh */
