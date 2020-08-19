//
// Namespace for collecting tools used in MC truth evaluation
// Original author: Dave Brown (LBNL) 8/10/2016
//
#ifndef TrkDiag_InfoMCStructHelper_hh
#define TrkDiag_InfoMCStructHelper_hh
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "RecoDataProducts/inc/KalSeed.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/GenInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
#include "TrkDiag/inc/CaloClusterInfoMC.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "MCDataProducts/inc/PrimaryParticle.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include <vector>
#include <functional>
namespace mu2e {

  class InfoMCStructHelper {

  private:
    art::InputTag _spctag;
    art::Handle<SimParticleCollection> _spcH;
    SimParticleTimeOffset _toff;
    double _mingood;

    void fillPriInfo(const PrimaryParticle& primary, GenInfo& priinfo);
    void fillGenInfo(const art::Ptr<GenParticle>& gp, GenInfo& geninfo);

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      
      fhicl::Atom<art::InputTag> spctag{Name("SimParticleCollectionTag"), Comment("InputTag for the SimParticleCollection"), art::InputTag()};
      fhicl::Sequence<art::InputTag> toff{Name("TimeMaps"), Comment("List of SimParticle time maps to use")};
      fhicl::Atom<double> mingood{Name("MinGoodMomFraction"), Comment("Minimum fraction of the true particle's momentum for a digi to be described as \"good\"")};
    };

    InfoMCStructHelper(const Config& conf) :
      _spctag(conf.spctag()), _toff(conf.toff()), _mingood(conf.mingood()) {  };

    void updateEvent(const art::Event& event) {
      event.getByLabel(_spctag,_spcH);
      _toff.updateMap(event);
    }

    const SimParticleTimeOffset& getTimeMaps() const { return _toff; }
    void fillTrkInfoMC(const KalSeedMC& kseedmc, TrkInfoMC& trkinfomc);
    void fillTrkInfoMCDigis(const KalSeedMC& kseedmc, TrkInfoMC& trkinfomc);
    void fillHitInfoMC(const KalSeedMC& kseedmc, TrkStrawHitInfoMC& tshinfomc, const TrkStrawHitMC& tshmc);
    void fillGenAndPriInfo(const KalSeedMC& kseedmc, const PrimaryParticle& primary, GenInfo& priinfo, GenInfo& geninfo);
    void fillTrkInfoMCStep(const KalSeedMC& kseedmc, TrkInfoMCStep& trkinfomcstep, std::vector<int> const& vids, double target_time);

    void fillHitInfoMCs(const KalSeedMC& kseedmc, std::vector<TrkStrawHitInfoMC>& tshinfomcs);
    void fillCaloClusterInfoMC(CaloClusterMC const& ccmc, CaloClusterInfoMC& ccimc);
  };
}

#endif
