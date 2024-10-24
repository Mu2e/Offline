// Ed Callaghan
// Decrypt a BigNumber using the Goldwasser-Micali cryposystem
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

// gmp
#include <gmp.h>
#include <gmpxx.h>

// mu2e
#include "Offline/Blinding/inc/BigNumber.hh"
#include "Offline/Blinding/inc/GMPRoutines.hh"

namespace mu2e{
  class GoldwasserMicaliDecrypter: public art::EDProducer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> encrypted_tag{
          fhicl::Name("BigNumber"),
          fhicl::Comment("art::InputTag of BigNumber to encrypt")
        };
        fhicl::Atom<std::string> prime_p{
          fhicl::Name("prime_p"),
          fhicl::Comment("Private key first prime (p)")
        };
        fhicl::Atom<std::string> prime_q{
          fhicl::Name("prime_q"),
          fhicl::Comment("Public key second prime (q)")
        };
      };

      using Parameters = art::EDProducer::Table<Config>;
      GoldwasserMicaliDecrypter(const Parameters&);
     ~GoldwasserMicaliDecrypter();
    protected:
      art::InputTag _encrypted_tag;
      BigNumber _prime_p;
      BigNumber _prime_q;
      mpz_t _p;
      mpz_t _q;
      mpz_t _N;
    private:
      void produce(art::Event&);
  };

  GoldwasserMicaliDecrypter::GoldwasserMicaliDecrypter(const Parameters& cfg):
      art::EDProducer(cfg),
      _encrypted_tag(cfg().encrypted_tag()),
      _prime_p(cfg().prime_p()),
      _prime_q(cfg().prime_q()){
    // allocate and initialize gmp representations
    mpz_inits(_p, _q, _N, NULL);
    mpz_set_str(_p, _prime_p.Buffer(), 10);
    mpz_set_str(_q, _prime_q.Buffer(), 10);
    mpz_mul(_N, _p, _q);

    // framework hooks
    this->consumes<BigNumber>(_encrypted_tag);
    this->produces<BigNumber>();
  }

  GoldwasserMicaliDecrypter::~GoldwasserMicaliDecrypter(){
    // deallocate gmp representation
    mpz_clears(_p, _q, _N, NULL);
  }

  void GoldwasserMicaliDecrypter::produce(art::Event& event){
    // locate unencrypted tag
    auto encrypted = event.getValidHandle<BigNumber>(_encrypted_tag);

    // local gmp initializations
    mpz_t gmp_decrypted;
    mpz_t gmp_encrypted;
    mpz_inits(gmp_decrypted, gmp_encrypted, NULL);
    mpz_set_str(gmp_encrypted, encrypted->Buffer(), 10);

    // apply the actual decryption scheme
    unsigned int set = gmp::gm_decrypt(gmp_encrypted, _p, _q, _N);

    // wrap encrypted tag as BigNumber, and store
    auto decrypted = std::make_unique<BigNumber>(std::to_string(set));
    event.put(std::move(decrypted));

    // lcoal gmp deallocations
    mpz_clears(gmp_decrypted, gmp_encrypted, NULL);
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GoldwasserMicaliDecrypter);
