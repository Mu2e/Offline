//
// Create STMHits from STMDigis
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
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

#include <utility>
// root
#include "TH1F.h"
#include "TTree.h"

#include "Offline/RecoDataProducts/inc/STMDigi.hh"
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
      fhicl::Atom<art::InputTag> stmDigisTag{ Name("stmDigisTag"), Comment("InputTag for STMDigiCollection")};
      //       fhicl::Atom<art::InputTag> adcNorm{ Name("adcNorm"), Comment("Normalization for the adc values")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit MakeSTMHits(const Parameters& conf);

  private:
    void produce(art::Event& e) override;

    art::InputTag _stmDigisTag;
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
    //    art::InputTag _adcNorm;
  };

  MakeSTMHits::MakeSTMHits(const Parameters& config )  :
    art::EDProducer{config},
    _stmDigisTag(config().stmDigisTag()),
    _stmEnergyCalib_h()
    //    _adcNorm(config().adcNorm())
    {
      consumes<STMDigiCollection>(_stmDigisTag);
      produces<STMHitCollection>();
    }

    void MakeSTMHits::produce(art::Event& event) {
    // create output
    unique_ptr<STMHitCollection> outputSTMHits(new STMHitCollection);
    auto digisHandle = event.getValidHandle<STMDigiCollection>(_stmDigisTag);

    if (digisHandle->size() > 0) {
      STMEnergyCalib const& stmEnergyCalib = _stmEnergyCalib_h.get(event.id()); // get calibration

      // Peek at first digi to get the channel since all digis in the same collection should be from the same detector
      STMChannel ch = digisHandle->at(0).channel();
      const auto& pars = stmEnergyCalib.calib(ch);
      //      std::cout << ch.name() << ": p0 = " << pars.p0 << ", p1 = " << pars.p1 << ", p2 = " << pars.p2 << std::endl;

      for (const auto& digi : *digisHandle) {
        int tdc = digi.trigTime();
        const std::vector<short int>& adc = digi.adcs();
        float time = tdc;
        auto uncalib_energy = *std::max_element(adc.begin(), adc.end());
        float energy = pars.p0 + pars.p1*uncalib_energy + pars.p2*uncalib_energy*uncalib_energy;
        STMHit stm_hit(time,energy);
        outputSTMHits->push_back(stm_hit);
      }
    }

    event.put(std::move(outputSTMHits));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMHits)
