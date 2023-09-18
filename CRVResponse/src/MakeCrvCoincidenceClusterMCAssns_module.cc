//
// Module that creates an Assns between a CrvCoincidenceCluster
// and its associated CrvCoincidenceClusterMC
//
// Can take many collections of CrvCoincidenceClusters
//
// Original author: Andy Edmonds 2023-05-02
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMCAssns.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

namespace mu2e {

  class MakeCrvCoincicdenceClusterMCAssns : public art::EDProducer {

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> crvCoincidenceTags{Name("crvCoincidenceTags"), Comment("art::InputTags for CrvCoincidenceClusterCollections")};
      fhicl::Sequence<art::InputTag> crvCoincidenceMCTags{Name("crvCoincidenceMCTags"), Comment("art::InputTags for CrvCoincidenceClusterMCCollections")};
    };
    typedef art::EDProducer::Table<Config> Parameters;

    explicit MakeCrvCoincicdenceClusterMCAssns(const Parameters& conf);
    virtual ~MakeCrvCoincicdenceClusterMCAssns() { }

    void beginJob();
    void produce(art::Event& e);

  private:

    // fhicl parameters
    std::vector<art::InputTag> _crvCoincidenceTags;
    std::vector<art::InputTag> _crvCoincidenceMCTags;
  };

  MakeCrvCoincicdenceClusterMCAssns::MakeCrvCoincicdenceClusterMCAssns(const Parameters& conf):
    art::EDProducer(conf),
    _crvCoincidenceTags(conf().crvCoincidenceTags()),
    _crvCoincidenceMCTags(conf().crvCoincidenceMCTags())
  {
    produces<CrvCoincidenceClusterMCAssns>();
  }

  void MakeCrvCoincicdenceClusterMCAssns::beginJob( ){  }

  void MakeCrvCoincicdenceClusterMCAssns::produce(art::Event& event) {

    auto outputAssns = std::make_unique<CrvCoincidenceClusterMCAssns>();
    for (size_t i_collection = 0; i_collection < _crvCoincidenceTags.size(); ++i_collection) {
      const auto& crvCoincidenceHandle = event.getValidHandle<CrvCoincidenceClusterCollection>(_crvCoincidenceTags.at(i_collection));
      const auto& crvCoincidenceMCHandle = event.getValidHandle<CrvCoincidenceClusterMCCollection>(_crvCoincidenceMCTags.at(i_collection));

      if (crvCoincidenceHandle->size() != crvCoincidenceMCHandle->size()) {
        throw cet::exception("MakeCrvCoincicdenceClusterMCAssns") << "CrvCoincidenceCollection and CrvCoincidenceMCCollection have different sizes (" << crvCoincidenceHandle->size() << " != " << crvCoincidenceMCHandle->size() << ")" << std::endl;
      }

      for(size_t i_crvCoinc = 0; i_crvCoinc != crvCoincidenceHandle->size(); ++i_crvCoinc) {
        auto crvCoincPtr = art::Ptr<CrvCoincidenceCluster>(crvCoincidenceHandle, i_crvCoinc);
        auto crvCoincMCPtr = art::Ptr<CrvCoincidenceClusterMC>(crvCoincidenceMCHandle, i_crvCoinc);
        outputAssns->addSingle(crvCoincPtr, crvCoincMCPtr);
      }
    }
    event.put(std::move(outputAssns));
  }
}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::MakeCrvCoincicdenceClusterMCAssns;
DEFINE_ART_MODULE(MakeCrvCoincicdenceClusterMCAssns)
