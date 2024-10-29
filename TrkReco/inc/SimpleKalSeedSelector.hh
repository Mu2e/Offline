//
//  Simple selector based on momentum and fit quality
//
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "Offline/TrkReco/inc/KalSeedSelector.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "fhiclcpp/types/Atom.h"

namespace mu2e {
  class SimpleKalSeedSelector : public KalSeedSelector {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<float> minmom{Name("MinMomentum"), Comment("Minimum fit momentum ")};
      fhicl::Atom<float> maxmom{Name("MaxMomentum"), Comment("Maximum fit momentum ")};
      fhicl::Atom<float> minfcon{Name("MinFitConsistency"), Comment("Minimum fit consistency ")};
      fhicl::Atom<unsigned> minnactive{Name("MinActiveHits"), Comment("Minimum # of active hits ")};
      fhicl::Atom<float> minsignhit{Name("MinDeltaNHitFraction"), Comment("Minimum difference in the fractional number of hits to consider significant")};
    };
    typedef art::ToolConfigTable<Config> Parameters;
    explicit SimpleKalSeedSelector(Parameters const& conf) :
      minmom_(conf().minmom()),
      maxmom_(conf().maxmom()),
      minfcon_(conf().minfcon()),
      minnactive_(conf().minnactive()),
      minsignhit_(conf().minsignhit())
    {}

    bool select(KalSeed const& kseed) const override;
    bool isBetter(KalSeed const& current,KalSeed const& test) const override;

  private:
    float minmom_, maxmom_;
    float minfcon_;
    unsigned minnactive_;
    float minsignhit_;
  };
}
DEFINE_ART_CLASS_TOOL(mu2e::SimpleKalSeedSelector)
