//
//  Utility class to print SimParticle
// 
#ifndef Print_inc_SimParticlePrinter_hh
#define Print_inc_SimParticlePrinter_hh

#include <cstring>
#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Print/inc/ProductPrinter.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class SimParticlePrinter : public ProductPrinter {
  public:

    struct Config : public ProductPrinter::Config {
      fhicl::Atom<float> pCut{ fhicl::Name("pCut"), 
	  fhicl::Comment("momentum cut on all particles"), -1 };
      fhicl::Atom<float> emPCut{ fhicl::Name("emPcut"), 
	  fhicl::Comment("momentum cut on EM particles"), -1 };
      fhicl::Atom<bool> primaryOnly{ fhicl::Name("primaryOnly"), 
	  fhicl::Comment("show only the primary particles"), false };
    };

    SimParticlePrinter():_pCut(-1),_emPCut(-1),_primaryOnly(false) {}
    SimParticlePrinter(const Config& conf):ProductPrinter(conf) { 
      _pCut = conf.pCut();
      _emPCut = conf.emPCut();
      _primaryOnly = conf.primaryOnly();
    }

    // do not print if p is below this cut
    void setPCut(double p) { _pCut = p; }
    // do not print e and gamma if p is below this cut
    void setEmPCut(double p) { _emPCut = p; }
    // print only particles marked primary
    void setPrimaryOnly(bool q) { _primaryOnly = q; }

    // all the ways to request a printout
    void Print(art::Event const& event,
	       std::ostream& os = std::cout) override;
    void Print(const art::Handle<SimParticleCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const art::ValidHandle<SimParticleCollection>& handle, 
	       std::ostream& os = std::cout);
    void Print(const SimParticleCollection& coll, 
	       std::ostream& os = std::cout);
    void Print(const art::Ptr<SimParticle>& ptr, 
	  int ind = -1, std::ostream& os = std::cout);
    void Print(const mu2e::SimParticle& obj, 
	  int ind = -1, std::size_t key = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:
    double _pCut;
    double _emPCut;
    bool _primaryOnly;

  };

}
#endif
