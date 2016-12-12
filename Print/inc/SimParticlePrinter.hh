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

    typedef std::vector<std::string> vecstr;

    SimParticlePrinter() { set( fhicl::ParameterSet() ); }
    SimParticlePrinter(const fhicl::ParameterSet& pset) { set(pset); }

    // tags to select which product instances to process
    void setTags(const vecstr& tags) { _tags = tags; }
    // usually customized in the subclass for items relevant to that product

    // methods to setup parameters
    // do not print if p is below this cut
    void setPCut(double p) { _pCut = p; }
    // do not print e and gamma if p is below this cut
    void setEmPCut(double p) { _emPCut = p; }
    // print only particles marked primary
    void setPrimaryOnly(bool q) { _primaryOnly = q; }

    // pset should contain a table called SimParticlePrinter
    void set(const fhicl::ParameterSet& pset);

    // the vector<string> list of inputTags
    const vecstr& tags() const {return _tags; }

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
	       int ind = -1, std::ostream& os = std::cout);

    void PrintHeader(const std::string& tag, 
		     std::ostream& os = std::cout);
    void PrintListHeader(std::ostream& os = std::cout);

  private:
    double _pCut;
    double _emPCut;
    bool _primaryOnly;
    vecstr _tags;

  };

}
#endif
