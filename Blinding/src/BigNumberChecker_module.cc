// Ed Callaghan
// Compare two BigNumbers for equality; useful for validating en-/de-cryption
// August 2024

// stl
#include <string>

// art
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"

// cetlib_except
#include "cetlib_except/exception.h"

// fhiclcpp
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

// mu2e
#include "Offline/Blinding/inc/BigNumber.hh"
#include "Offline/Blinding/inc/GMPRoutines.hh"

namespace mu2e{
  class BigNumberChecker: public art::EDAnalyzer{
    public:
      struct Config{
        fhicl::Atom<art::InputTag> lhs_tag{
          fhicl::Name("BigNumber_a"),
          fhicl::Comment("art::InputTag of BigNumber to compare to")
        };
        fhicl::Atom<art::InputTag> rhs_tag{
          fhicl::Name("BigNumber_b"),
          fhicl::Comment("art::InputTag of BigNumber to compare")
        };
      };

      using Parameters = art::EDAnalyzer::Table<Config>;
      BigNumberChecker(const Parameters&);
    protected:
      art::InputTag _lhs_tag;
      art::InputTag _rhs_tag;
    private:
      void analyze(const art::Event&);
  };

  BigNumberChecker::BigNumberChecker(const Parameters& cfg):
      art::EDAnalyzer(cfg),
      _lhs_tag(cfg().lhs_tag()),
      _rhs_tag(cfg().rhs_tag()){
    /**/
  }

  void BigNumberChecker::analyze(const art::Event& event){
    // locate BigNumbers
    auto lhs_handle = event.getValidHandle<BigNumber>(_lhs_tag);
    auto rhs_handle = event.getValidHandle<BigNumber>(_rhs_tag);

    auto lhs = *lhs_handle;
    auto rhs = *rhs_handle;
    if (lhs != rhs){
      std::string msg = "BigNumber mismatch: "
                      + lhs.String() + " != " + rhs.String();
      throw cet::exception("BigNumberChecker") << msg << std::endl;
    }
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::BigNumberChecker);
