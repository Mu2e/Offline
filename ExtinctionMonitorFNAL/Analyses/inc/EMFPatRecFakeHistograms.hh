// Track finding efficiency histograms.
//
// Andrei Gaponenko
//

#ifndef ExtinctionMonitorFNAL_Analyses_EMFPatRecFakeHistograms_hh
#define ExtinctionMonitorFNAL_Analyses_EMFPatRecFakeHistograms_hh

#include <string>

#include "art/Framework/Core/FindMany.h"
#include "art/Persistency/Common/Ptr.h"

#include "boost/noncopyable.hpp"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParamCollection.hh"
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

    class EMFPatRecFakeHistograms :  private boost::noncopyable {
    public:

      explicit EMFPatRecFakeHistograms(unsigned cutMinCommonClusters);

      // Book histograms in the subdirectory, given by the relativePath; that path is
      // relative to the root TFileDirectory for the current module.
      // The default is to use the current module's TFileDirectory
      void book(const ExtMon& extmon, const std::string& relativePath="");

      // Book histograms in the specified TFileDirectory.
      void book(const ExtMon& extmon, art::TFileDirectory& tfdir);

      // A helper class to use per-event FinMany object
      // instead of creating one per fill call.
      class Fillable {
        EMFPatRecFakeHistograms *parent_;
        art::FindMany<SimParticle,ExtMonFNALTrkMatchInfo> particleFinder_;
        unsigned multiplicity_; // per-event variable
        std::set<const SimParticle*> usedParticles_; // already assigned to tracks
      public:

        Fillable(EMFPatRecFakeHistograms *parent,
                 const art::Handle<ExtMonFNALTrkParamCollection>& tracks,
                 const art::Event& event,
                 const art::InputTag& trkTruthTag,
                 unsigned multiplicity)
          : parent_(parent)
          , particleFinder_(tracks, event, trkTruthTag)
          , multiplicity_(multiplicity)
        {}

        void fill(unsigned denominatorTrackIndex);
      };

      Fillable fillable(const art::Handle<ExtMonFNALTrkParamCollection>& tracks,
                        const art::Event& event,
                        const art::InputTag& trkTruthTag,
                        unsigned multiplicity)  {
        return Fillable(this, tracks, event, trkTruthTag, multiplicity);
      }

    private:
      friend class Fillable;
      unsigned cutMinCommonClusters_; // inclusive

      ExtMon *extmon_; // non-owning

      TEfficiency* fracDups_;
      TEfficiency* fracInteractions_;
      TEfficiency* fracFakes_;
    };

  } // end namespace ExtMonFNAL
} // end namespace mu2e

#endif /* ExtinctionMonitorFNAL_Analyses_EMFPatRecFakeHistograms_hh */
