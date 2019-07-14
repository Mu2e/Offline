//Hasung Song (2018)

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib_except/exception.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"

#include "TF1.h"

namespace mu2e {

  class RMCPhys {
    
    private :
      art::InputTag input_;
      double kmax_;
      int internalconversion_;
      int verbosityLevel_;
    public :
    RMCPhys(const fhicl::ParameterSet& pset)
      : input_(pset.get<std::string>("inputModule","compressDigiMCs") )
      , kmax_(pset.get<double>("kinematic_endpoint") )
      , internalconversion_(pset.get<int>("internalConversion",0))
      , verbosityLevel_(pset.get<int>("verbosityLevel", 0 ) ) { }
    double x;
    TF1 *RMCSpectrum = new TF1("RMC","(1-2*x+2*(x*x))*x*(1-x)*(1-x)",0,1);
    double norm = 1./RMCSpectrum->Integral(57./kmax_,1);
    double weight(const art::Event& evt) {
      auto genColl = evt.getValidHandle<GenParticleCollection>( input_ );
      double energy = 0;
      double wt = 0;
      if (internalconversion_ > 0) {
        for ( const auto& i: *genColl ) {
          if (i.generatorId().id() ==  42) { //22-internalRPC, 41-RMC, 42-internalRMC
            energy += i.momentum().e();
            if (verbosityLevel_ ==1) {
              std::cout << "Particle id:" << i.pdgId() << "& momentum: " << i.momentum().e() << std::endl;
            }
          }
        } 
      }
      else {
        for ( const auto& i: *genColl ) {
          if (i.pdgId() == 22 && i.generatorId().id() == 41 ) {
            energy = i.momentum().e(); 
          }
        }
      }
      if (energy > kmax_) wt = 0;
      else {
        x = energy/kmax_;
        wt = norm*(1-2*x+2*(x*x))*x*(1-x)*(1-x) ; //RMCSpectrum1;
      }
      if (verbosityLevel_ == 1) {
        std::cout << "Norm: " << norm << std::endl;
        std::cout << "Energy: " << energy << std::endl;
        std::cout << "Weight:" << wt << std::endl;
      
      }
      return wt;
    };
  };
}
