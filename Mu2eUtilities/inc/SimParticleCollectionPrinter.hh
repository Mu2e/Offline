// Andrei Gaponenko, 2013

#ifndef Mu2eUtilities_inc_SimParticleCollectionPrinter_hh
#define Mu2eUtilities_inc_SimParticleCollectionPrinter_hh

#include <ostream>
#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {
  class SimParticleCollectionPrinter {
    std::string prefix_;
    bool enabled_;
    bool primariesOnly_;
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string> prefix {Name("prefix"),
          Comment("A string to print at the start of each new line, before a SimParticle printout"),
          ""
          };

      fhicl::Atom<bool> primariesOnly {Name("primariesOnly"),
          Comment("An option to print just primary particles; the default is all particles."),
          false};

      fhicl::Atom<bool> enabled {Name("enabled"), true};
    };


    explicit SimParticleCollectionPrinter(const Config& conf);

    // A constructor to create an instance that is disabled by default.
    // For use with fhicl::OptionalTable.
    SimParticleCollectionPrinter();

    static std::ostream& print(std::ostream& os, const SimParticle& particle);

    std::ostream& print(std::ostream& os, const SimParticleCollection& particles) const;

    // This can be used in module constructors to initialize
    // printer instances that are not enabled by default.
    static fhicl::ParameterSet defaultPSet();
  };

}

#endif/*Mu2eUtilities_inc_SimParticleCollectionPrinter_hh*/
