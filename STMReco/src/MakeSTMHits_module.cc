//
// Create STMHits from STMPHDigis
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

#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
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
      fhicl::Atom<art::InputTag> stmPHDigisTag{ Name("stmPHDigisTag"), Comment("InputTag for STMPHDigiCollection")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit MakeSTMHits(const Parameters& conf);

  private:
    void produce(art::Event& e) override;

    art::ProductToken<STMPHDigiCollection> _stmPHDigisToken;
    STMChannel _channel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
  };

  MakeSTMHits::MakeSTMHits(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmPHDigisToken(consumes<STMPHDigiCollection>(config().stmPHDigisTag()))
    ,_channel(STMChannel::LaBr)
    ,_stmEnergyCalib_h()
    {
      produces<STMHitCollection>();
    }

    void MakeSTMHits::produce(art::Event& event) {
    // create output
    unique_ptr<STMHitCollection> outputSTMHits(new STMHitCollection);
    auto phDigisHandle = event.getValidHandle(_stmPHDigisToken);

    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get calibration

    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);
    const auto& pars = stmEnergyCalib.calib(_channel);

    for (const auto& ph_digi : *phDigisHandle) {
      auto uncalib_time = ph_digi.time();
      auto uncalib_energy = ph_digi.energy();
      float time = uncalib_time*nsPerCt;
      float energy = pars.p0 + pars.p1*uncalib_energy + pars.p2*uncalib_energy*uncalib_energy;

      STMHit stm_hit(time,energy);
      outputSTMHits->push_back(stm_hit);
    }

    event.put(std::move(outputSTMHits));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMHits)
