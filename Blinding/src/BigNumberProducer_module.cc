// Ed Callaghan
// Produce constant BigNumbers
// August 2024

// stl
#include <string>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/Blinding/inc/BigNumber.hh"

namespace mu2e{
  class BigNumberProducer: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<std::string> value{
          fhicl::Name("value"),
          fhicl::Comment("Value to embed")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      BigNumberProducer(const Parameters&);
    protected:
      std::string _digits;
    private:
      void produce(art::Event&);
  };

  BigNumberProducer::BigNumberProducer(const Parameters& config):
      art::EDProducer(config),
      _digits(config().value()){
    this->produces<BigNumber>();
  }

  void BigNumberProducer::produce(art::Event& event){
    auto product = std::make_unique<BigNumber>(_digits);
    event.put(std::move(product));
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::BigNumberProducer);
