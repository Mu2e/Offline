// Convert proton bunch time marker from data to MC data product
// To be used for mixing simulated signal into real data from detector
//
// Ed Callaghan, 2024

// art
#include <art/Framework/Core/EDProducer.h>
#include <art/Framework/Principal/Event.h>

// canvas
#include <canvas/Utilities/InputTag.h>

// fhiclcpp
#include <fhiclcpp/types/Atom.h>

// mu2e
#include <Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh>
#include <Offline/RecoDataProducts/inc/ProtonBunchTime.hh>

namespace mu2e{
  class ProtonBunchTimeMCFromProtonBunchTime: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> pbttag{
          fhicl::Name("PBTTag"),
          fhicl::Comment("ProtonBunchTime provenance"),
          "" // default, eventually to specify "from data"
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      explicit ProtonBunchTimeMCFromProtonBunchTime(Parameters const& config);
    protected:
      /**/
    private:
      void produce(art::Event& event) override;
      art::InputTag pbttag_;
  };

  ProtonBunchTimeMCFromProtonBunchTime::ProtonBunchTimeMCFromProtonBunchTime(Parameters const& config): EDProducer(config), pbttag_(config().pbttag()){
    this->produces<ProtonBunchTimeMC>();
  }

  void ProtonBunchTimeMCFromProtonBunchTime::produce(art::Event& event){
    auto pbt = event.getHandle<ProtonBunchTime>(pbttag_);
    std::unique_ptr<ProtonBunchTimeMC> pbtmc(new ProtonBunchTimeMC);
    pbtmc->pbtime_ = pbt->pbtime_;
    event.put(std::move(pbtmc));
  }
}

DEFINE_ART_MODULE(mu2e::ProtonBunchTimeMCFromProtonBunchTime);
