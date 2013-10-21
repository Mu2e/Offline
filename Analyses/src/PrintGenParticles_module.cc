// Prints out all GenParticles in a collection.
//
// $Id: PrintGenParticles_module.cc,v 1.2 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Andrei Gaponenko
//

#include <iostream>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

namespace mu2e {

  class PrintGenParticles : public art::EDAnalyzer {
  public:

    explicit PrintGenParticles(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e);

  private:

    // The two strings specify what collection to print
    std::string moduleLabel_;
    std::string instanceName_;
  };

  PrintGenParticles::PrintGenParticles(const fhicl::ParameterSet& pset) : 
    art::EDAnalyzer(pset),
    moduleLabel_(pset.get<std::string>("inputModuleLabel")),
    instanceName_(pset.get<std::string>("inputInstanceName"))
  {}

  void PrintGenParticles::analyze(const art::Event& event) {

    art::Handle<GenParticleCollection> ih;
    event.getByLabel(moduleLabel_, instanceName_, ih);
    const GenParticleCollection& gp(*ih);

    std::cout<<"PrintGenParticles: begin printing collection moduleLabel = \""<<moduleLabel_
             <<"\", instanceName = \""<<instanceName_<<"\""
             <<std::endl;
    for (GenParticleCollection::const_iterator i = gp.begin(); i != gp.end(); ++i) {
      std::cout<<"PrintGenParticles: "<<*i<< std::endl;
    }
    std::cout<<"PrintGenParticles:   end printing collection moduleLabel = \""<<moduleLabel_
             <<"\", instanceName = \""<<instanceName_<<"\""
             <<std::endl;

  } // analyze()

}  // end namespace mu2e

// Register the module with the framework
//using mu2e::PrintGenParticles;
DEFINE_ART_MODULE(mu2e::PrintGenParticles);
