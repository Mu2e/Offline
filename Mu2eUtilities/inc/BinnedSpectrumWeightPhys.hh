// A. Edmonds based on RMCPhys.hh by Hasung Song (2018)

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
#include "fhiclcpp/ParameterSet.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"

namespace mu2e {

  class BinnedSpectrumWeightPhys {
    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERGY, MOMENTUM };
    
  private :
    fhicl::ParameterSet psphys_;
    art::InputTag input_;
    int verbosityLevel_;
    PDGCode::type genPdg_;
    GenId genId_;
    SpectrumVar       spectrumVariable_;
    BinnedSpectrum spectrum_;

    static SpectrumVar parseSpectrumVar(const std::string& name) {
      if (name == "totalEnergy"  )  return TOTAL_ENERGY;
      if (name == "kineticEnergy")  return KINETIC_ENERGY;
      if (name == "momentum"     )  return MOMENTUM;
      throw cet::exception("BADCONFIG")<<"BinnedSpectrumWeightPhys: unknown spectrum variable "<<name<<"\n";
    }

  public :

    BinnedSpectrumWeightPhys(const fhicl::ParameterSet& pset)
      : psphys_(pset.get<fhicl::ParameterSet>("physics"))
      , input_(pset.get<std::string>("genParticleTag","compressDigiMCs") )
      , verbosityLevel_(pset.get<int>("verbosityLevel", 0 ) )
      , genPdg_(PDGCode::type(pset.get<int>("genParticlePdgId",0)) )
      , genId_(GenId::findByName(pset.get<std::string>("genParticleGenId")))
      , spectrumVariable_(parseSpectrumVar(psphys_.get<std::string>("spectrumVariable")))
      , spectrum_(BinnedSpectrum(psphys_))

    {
    }


    double weight(const art::Event& evt) {
      auto genColl = evt.getValidHandle<GenParticleCollection>( input_ );
      double sampleVal = 0;
      if (spectrumVariable_ == TOTAL_ENERGY){
        for ( const auto& i: *genColl ) {
          if ((i.pdgId() == genPdg_ || genPdg_ == PDGCode::null) && i.generatorId() == genId_ ) {
            sampleVal += i.momentum().e(); 
          }
        }
      }else if (spectrumVariable_ == KINETIC_ENERGY){
        for ( const auto& i: *genColl ) {
          if ((i.pdgId() == genPdg_ || genPdg_ == PDGCode::null) && i.generatorId() == genId_ ) {
            sampleVal += i.momentum().e() - i.momentum().restMass(); 
          }
        }
      }else{
        CLHEP::Hep3Vector mom(0,0,0);
        for ( const auto& i: *genColl ) {
          if ((i.pdgId() == genPdg_ || genPdg_ == PDGCode::null) && i.generatorId() == genId_ ) {
            mom += i.momentum().vect(); 
          }
        }
        sampleVal = mom.mag();
      }

      if (sampleVal > spectrum_.getXMax() || sampleVal < spectrum_.getXMin())
        return 0;

      size_t i_bin = (sampleVal - spectrum_.getXMin()) / spectrum_.getBinWidth();
      return spectrum_.getPDF(i_bin);
    };
  };
}
