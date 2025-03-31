// Ed Callaghan
// Encrypt a BigNumber using the Goldwasser-Micali cryposystem
// August 2024

// stl
#include <string>

// art
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// gmp
#include <gmp.h>
#include <gmpxx.h>

// mu2e
#include "Offline/Blinding/inc/BigNumber.hh"
#include "Offline/Blinding/inc/GMPRoutines.hh"
#include "Offline/SeedService/inc/SeedService.hh"

namespace mu2e{
  class GoldwasserMicaliEncrypter: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> unencrypted_tag{
          fhicl::Name("BigNumber"),
          fhicl::Comment("art::InputTag of BigNumber to encrypt")
        };
        fhicl::Atom<std::string> modulus{
          fhicl::Name("modulus"),
          fhicl::Comment("Public key modulus (N)")
        };
        fhicl::Atom<std::string> nonresidue{
          fhicl::Name("nonresidue"),
          fhicl::Comment("Public key nonresidue (x)")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      GoldwasserMicaliEncrypter(const Parameters&);
     ~GoldwasserMicaliEncrypter();
    protected:
      art::InputTag _unencrypted_tag;
      BigNumber _modulus;
      BigNumber _nonresidue;
      gmp_randstate_t _gmp_rand_state;
      mpz_t _x;
      mpz_t _N;
    private:
      void produce(art::Event&);
  };

  GoldwasserMicaliEncrypter::GoldwasserMicaliEncrypter(const Parameters& cfg):
      art::EDProducer(cfg),
      _unencrypted_tag(cfg().unencrypted_tag()),
      _modulus(cfg().modulus()),
      _nonresidue(cfg().nonresidue()){
    // allocate and initialize gmp representations
    mpz_inits(_x, _N, NULL);
    mpz_set_str(_x, _nonresidue.Buffer(), 10);
    mpz_set_str(_N, _modulus.Buffer(), 10);

    // seed gmp random state from on high
    auto service = art::ServiceHandle<SeedService>();
    auto seed = static_cast<unsigned int>(service->getSeed());
    gmp_randinit_default(_gmp_rand_state);
    gmp_randseed_ui(_gmp_rand_state, seed);

    // framework hooks
    this->consumes<BigNumber>(_unencrypted_tag);
    this->produces<BigNumber>();
  }

  GoldwasserMicaliEncrypter::~GoldwasserMicaliEncrypter(){
    // deallocate gmp representation
    mpz_clears(_x, _N, NULL);
  }

  void GoldwasserMicaliEncrypter::produce(art::Event& event){
    // locate unencrypted tag
    auto unencrypted = event.getValidHandle<BigNumber>(_unencrypted_tag);

    // local gmp initializations
    mpz_t gmp_unencrypted;
    mpz_t gmp_encrypted;
    mpz_inits(gmp_unencrypted, gmp_encrypted, NULL);
    mpz_set_str(gmp_unencrypted, unencrypted->Buffer(), 10);

    // apply the actual encryption scheme
    unsigned int set = 0;
    if (!unencrypted->IsZero()){
      set = 1;
    }
    gmp::gm_encrypt(gmp_encrypted, set, _x, _N, _gmp_rand_state);

    // wrap encrypted tag as BigNumber, and store
    auto encrypted = std::make_unique<BigNumber>(mpz_class(gmp_encrypted).get_str());
    event.put(std::move(encrypted));

    // lcoal gmp deallocations
    mpz_clears(gmp_unencrypted, gmp_encrypted, NULL);
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GoldwasserMicaliEncrypter);
