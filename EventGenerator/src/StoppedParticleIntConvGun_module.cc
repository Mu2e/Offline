// Andrei Gaponenko, 2013

// C++ includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

// cetlib includes
#include "cetlib/exception.h"

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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/Random2Dpair.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"

namespace mu2e {

  //================================================================
  class StoppedParticleIntConvGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    double elow_; // BinnedSpectrum does not store emin and emax reliably
    double ehi_;
    BinnedSpectrum spectrum_;
    static BinnedSpectrum parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                             double *elow,
                                             double *ehi);

    int verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    RandomUnitSphere randomUnitSphere_;

    RootTreeSampler<IO::StoppedParticleTauNormF> stops_;

    bool doHistograms_;

    double generateEnergy();

    TH1F* _hmomentum;
    TH1F* _hElecMom ;
    TH1F* _hPosiMom ;
    
  public:
    explicit StoppedParticleIntConvGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedParticleIntConvGun::StoppedParticleIntConvGun(const fhicl::ParameterSet& pset)
    : psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , elow_()
    , ehi_()
    , spectrum_(parseSpectrumShape(psphys_, &elow_, &ehi_))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomUnitSphere_(eng_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("pionStops"))
    , doHistograms_( pset.get<bool>("doHistograms",true) )
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::EventWeight>();

    if(verbosityLevel_ > 0) {
      std::cout<<"StoppedParticleIntConvGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StoppedParticleIntConvGun: producing electron-positron pair " << std::endl;
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "StoppedParticleIntConvGun" );

      _hmomentum = tfdir.make<TH1F>( "hmomentum", "Produced photon momentum"  , 100,  40.,  140.  );
      _hElecMom  = tfdir.make<TH1F>( "hElecMom" , "Produced electron momentum", 140,  0. ,  140.  );
      _hPosiMom  = tfdir.make<TH1F>( "hPosiMom" , "Produced positron momentum", 140,  0. ,  140.  );

    }

  }

  //================================================================
  BinnedSpectrum
  StoppedParticleIntConvGun::parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                                 double *elow,
                                                 double *ehi)
  {
    BinnedSpectrum res;

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    if (spectrumShape == "Bistirlich") {
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      res.initialize<PionCaptureSpectrum>( *elow, *ehi, psphys.get<double>("spectrumResolution") );
    }
    else if (spectrumShape == "flat") {
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      res.initialize<SimpleSpectrum>(*elow, *ehi, *ehi-*elow, SimpleSpectrum::Spectrum::Flat );
    }
    else {
      throw cet::exception("BADCONFIG")
        << "StoppedParticlePionGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }

    return res;
  }

  //================================================================
  void StoppedParticleIntConvGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const double energy = generateEnergy();

    // Need mass of electron
    static const double massE = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus).ref().mass().value();

    // Uses energy above as photon energy, assuming distribution created by 
    static Random2Dpair< PionCaptureSpectrum > random2dPair( eng_, 2*massE, energy, -1., 1. );
        
    const auto xyPair          = random2dPair.fire( energy );
    const auto elecPosiVectors = PionCaptureSpectrum::getElecPosiVectors( energy, xyPair.first, xyPair.second ); 
       
    // Add particles to list
    output->emplace_back( PDGCode::e_minus, GenId::internalRPC, pos, elecPosiVectors.first , stop.t );
    output->emplace_back( PDGCode::e_plus , GenId::internalRPC, pos, elecPosiVectors.second, stop.t );

    event.put(std::move(output));

    // Calculate survival probability
    const double weight = exp(-stop.tauNormalized);
    std::unique_ptr<EventWeight> pw(new EventWeight(weight));
    event.put(std::move(pw));

    if ( doHistograms_ ) {
      _hmomentum->Fill( energy );
      _hElecMom ->Fill( elecPosiVectors.first .vect().mag() );
      _hPosiMom ->Fill( elecPosiVectors.second.vect().mag() );
    }

  }

  //================================================================
  double StoppedParticleIntConvGun::generateEnergy() {
    return elow_ + (ehi_ - elow_)*randSpectrum_.fire();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticleIntConvGun);
