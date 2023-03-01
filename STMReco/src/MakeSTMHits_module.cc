//
// Create STMHits from STMMWDDigis
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

#include <utility>
// root
#include "TH1F.h"
#include "TTree.h"

#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
#include "Offline/RecoDataProducts/inc/STMHit.hh"

// C++
#include <vector>

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class MakeSTMHits : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> stmMWDDigisTag{ Name("stmMWDDigisTag"), Comment("InputTag for STMMWDDigiCollection")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit MakeSTMHits(const Parameters& conf);

  private:
    void produce(art::Event& e) override;

    art::ProductToken<STMMWDDigiCollection> _stmMWDDigisToken;
    STMChannel _channel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
  };

  MakeSTMHits::MakeSTMHits(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmMWDDigisToken(consumes<STMMWDDigiCollection>(config().stmMWDDigisTag()))
    ,_channel(STMUtils::getChannel(config().stmMWDDigisTag()))
    ,_stmEnergyCalib_h()
    {
      produces<STMHitCollection>();
    }

    void MakeSTMHits::produce(art::Event& event) {
    // create output
    unique_ptr<STMHitCollection> outputSTMHits(new STMHitCollection);
    auto mwdDigisHandle = event.getValidHandle(_stmMWDDigisToken);

    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get calibration

    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);
    const auto& pars = stmEnergyCalib.calib(_channel);

    for (const auto& mwd_digi : *mwdDigisHandle) {
      auto uncalib_time = mwd_digi.time();
      auto uncalib_energy = mwd_digi.energy();
      float time = uncalib_time*nsPerCt;
      float energy = pars.p0 + pars.p1*uncalib_energy + pars.p2*uncalib_energy*uncalib_energy;

      STMHit stm_hit(time,energy);
      outputSTMHits->push_back(stm_hit);
    }

    event.put(std::move(outputSTMHits));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMHits)
