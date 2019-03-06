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
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"

namespace mu2e {

  class BinnedSpectrumWeightPhys {
    
  private :
    art::InputTag input_;
    mu2e::BinnedSpectrum spectrum_;
    int verbosityLevel_;
    double elow_;
    double ehi_;
    PDGCode::type genPdg_;
    GenId genId_;

  public :
    BinnedSpectrumWeightPhys(const fhicl::ParameterSet& pset)
      : input_(pset.get<std::string>("genParticleTag","compressDigiMCs") )
      , verbosityLevel_(pset.get<int>("verbosityLevel", 0 ) )
      , genPdg_(PDGCode::type(pset.get<int>("genParticlePdgId")) )
      , genId_(GenId::findByName(pset.get<std::string>("genParticleGenId")) ){

      std::string spectrumFileName = pset.get<std::string>("spectrumFileName");
      spectrum_.initialize(loadTable<2>( ConfigFileLookupPolicy()( spectrumFileName )));
      elow_ = spectrum_.getAbscissa(0);
      ehi_  = spectrum_.getAbscissa(spectrum_.getNbins()-1) + spectrum_.getBinWidth();
      if (pset.get<bool>("BinCenter", false)) {
        elow_ -= spectrum_.getBinWidth()/2;
        ehi_  -= spectrum_.getBinWidth()/2;
      }
      if(elow_ < 0.0) {
	throw cet::exception("BADCONFIG")
	  << "StoppedParticleReactionGun: negative energy endpoint "<< elow_ <<"\n";
      }
    }

    double weight(const art::Event& evt) {
      auto genColl = evt.getValidHandle<GenParticleCollection>( input_ );
      double energy = 0;
      double wt = 0;
      for ( const auto& i: *genColl ) {
        if (i.pdgId() == genPdg_ && i.generatorId() == genId_ ) {
          energy = i.momentum().e(); 
        }
      }

      if (energy > ehi_) {
	wt = 0;
      }
      else {
	size_t i_bin = (energy - elow_) / spectrum_.getBinWidth();
	wt = spectrum_.getPDF(i_bin);
      }
      return wt;
    };
  };
}
