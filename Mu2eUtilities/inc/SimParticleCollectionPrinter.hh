// Andrei Gaponenko, 2013

#ifndef Mu2eUtilities_inc_SimParticleCollectionPrinter_hh
#define Mu2eUtilities_inc_SimParticleCollectionPrinter_hh

#include <ostream>
#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {
  class SimParticleCollectionPrinter {
    std::string prefix_;
    bool enabled_;
    bool primariesOnly_;
  public:

    // Could configure printout format via pset
    explicit SimParticleCollectionPrinter(const fhicl::ParameterSet& pset);

    static std::ostream& print(std::ostream& os, const SimParticle& particle);

    std::ostream& print(std::ostream& os, const SimParticleCollection& particles) const;

    // This can be used in module constructors to initialize
    // printer instances that are not enabled by default.
    static fhicl::ParameterSet defaultPSet();
  };

}

#endif/*Mu2eUtilities_inc_SimParticleCollectionPrinter_hh*/
