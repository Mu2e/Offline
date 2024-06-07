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
      fhicl::Atom<double> minmom{Name("MinimumMomentum"), Comment("Minimum fit momentum ")};
      fhicl::Atom<double> maxmom{Name("MaximumMomentum"), Comment("Maximum fit momentum ")};
      fhicl::Atom<double> minfcon{Name("MinimumFitConsistency"), Comment("Minimum fit consistency ")};
      fhicl::Atom<double> minsignhit{Name("MinDeltaNHitFraction"), Comment("Minimum difference in the fractional number of hits to consider significant")};
    };
    typedef art::ToolConfigTable<Config> Parameters;
    explicit SimpleKalSeedSelector(Parameters const& conf) :
      minmom_(conf().minmom()),
      maxmom_(conf().maxmom()),
      minfcon_(conf().minfcon()),
      minsignhit_(conf().minsignhit())
    {}

    bool select(KalSeed const& kseed) const override;
    bool isBetter(KalSeed const& current,KalSeed const& test) const override;

  private:
    double minmom_, maxmom_;
    double minfcon_;
    double minsignhit_;
  };
}
DEFINE_ART_CLASS_TOOL(mu2e::SimpleKalSeedSelector)
