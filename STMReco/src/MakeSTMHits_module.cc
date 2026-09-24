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
      fhicl::Atom<art::InputTag> stmPHDigisMapTag{ Name("stmPHDigisMapTag"), Comment("InputTag for STMPHDigiCollectionMap")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit MakeSTMHits(const Parameters& conf);

  private:
    void produce(art::Event& e) override;

    art::ProductToken<STMPHDigiCollectionMap> _stmPHDigisMapToken;
    STMChannel _channel;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
  };

  MakeSTMHits::MakeSTMHits(const Parameters& config )  :
    art::EDProducer{config}
    ,_stmPHDigisMapToken(consumes<STMPHDigiCollectionMap>(config().stmPHDigisMapTag()))
    ,_channel(STMUtils::getChannel(config().stmPHDigisMapTag()))
    ,_stmEnergyCalib_h()
    {
      produces<STMHitCollectionMap>();
    }

    void MakeSTMHits::produce(art::Event& event) {
    // create output
    unique_ptr<STMHitCollectionMap> outputSTMHitsMap(new STMHitCollectionMap);
    auto phDigisHandle = event.getValidHandle(_stmPHDigisMapToken);

    STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get calibration
    const auto nsPerCt = stmEnergyCalib.nsPerCt(_channel);
    const auto& pars = stmEnergyCalib.calib(_channel);

    for (const auto& mu2e_evt : *phDigisHandle) {

      const auto& stm_evt_header = mu2e_evt.first; // get header information
      const auto& ph_digis = mu2e_evt.second;

      for (const auto& ph_digi : ph_digis) {
        auto uncalib_time = ph_digi.time();
        auto uncalib_energy = ph_digi.energy();
        // make a hit
        float time = uncalib_time*nsPerCt;
        float energy = pars.p0 + pars.p1*uncalib_energy + pars.p2*uncalib_energy*uncalib_energy;
        // store hit
        STMHit stm_hit(time,energy);
        // Add to map
        (*outputSTMHitsMap)[stm_evt_header].push_back(stm_hit);
      }
    }

    event.put(std::move(outputSTMHitsMap));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMHits)
