// Andrei Gaponenko, 2013

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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "DataProducts/inc/PDGCode.hh"
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
#include "TH1.h"
#include "TH2.h"

using namespace std::string_literals;


namespace mu2e {

  //================================================================
  class StoppedParticleIntConvGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    double elow_; // BinnedSpectrum does not store emin and emax reliably
    double ehi_;
    // Spectrum must be constructed after elow_ and ehi_.
    BinnedSpectrum spectrum_;
    int            verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;

    CLHEP::RandGeneral* randSpectrum_;
    RandomUnitSphere    randomUnitSphere_;
    CLHEP::RandFlat     randomFlat_;

    PionCaptureSpectrum pionCaptureSpectrum_;

    RootTreeSampler<IO::StoppedParticleTauNormF> stops_;

    bool   doHistograms_;

    TH1F* _hmomentum{nullptr};
    TH1F* _hElecMom {nullptr};
    TH1F* _hPosiMom {nullptr};
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;    		// M(ee)/E(gamma)
    TH1F* _hy;				// splitting function
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    double generateEnergy();
    void   parseSpectrumShape(const fhicl::ParameterSet& psphys);

  public:
    explicit StoppedParticleIntConvGun(const fhicl::ParameterSet& pset);
    ~StoppedParticleIntConvGun();
    void produce(art::Event& event) override;
  };

  //================================================================
  StoppedParticleIntConvGun::StoppedParticleIntConvGun(const fhicl::ParameterSet& pset):
    psphys_ (pset.get<fhicl::ParameterSet>("physics"))
    , elow_   (psphys_.get<double>("elow"))
    , ehi_    (psphys_.get<double>("ehi" ))
    , verbosityLevel_{pset.get<int>("verbosityLevel", 0)}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randomUnitSphere_(eng_)
    , randomFlat_(eng_)
    , pionCaptureSpectrum_(&randomFlat_,&randomUnitSphere_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("pionStops"))
    , doHistograms_{pset.get<bool>("doHistograms", true)}
  {
    produces<GenParticleCollection>();
    produces<EventWeight>();

    if (verbosityLevel_ > 0) {
      std::cout << "StoppedParticleIntConvGun: using = "
                << stops_.numRecords()
                << " stopped particles\n"
                << "StoppedParticleIntConvGun: producing electron-positron pair " << std::endl;
    }

    // initialize binned spectrum - this needs to be done right
    parseSpectrumShape(psphys_);
    randSpectrum_ = new CLHEP::RandGeneral(eng_, spectrum_.getPDF(), spectrum_.getNbins());

    if (doHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("StoppedParticleIntConvGun");

      _hmomentum = tfdir.make<TH1F>("hmomentum", "Produced photon momentum"  , 100,  40., 140.);
      _hElecMom  = tfdir.make<TH1F>("hElecMom" , "Produced electron momentum", 140,  0. , 140.);
      _hPosiMom  = tfdir.make<TH1F>("hPosiMom" , "Produced positron momentum", 140,  0. , 140.);
      _hMee      = tfdir.make<TH1F>("hMee"     , "M(e+e-) "           , 200,0.,200.);
      _hMeeVsE   = tfdir.make<TH2F>("hMeeVsE"  , "M(e+e-) vs E"       , 200,0.,200.,200,0,200);
      _hMeeOverE = tfdir.make<TH1F>("hMeeOverE", "M(e+e-)/E "         , 200, 0.,1);
      _hy        = tfdir.make<TH1F>("hy"       , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
    }
  }

  //-----------------------------------------------------------------------------
  StoppedParticleIntConvGun::~StoppedParticleIntConvGun() {
    delete randSpectrum_;
  }

//-----------------------------------------------------------------------------
  void StoppedParticleIntConvGun::parseSpectrumShape(const fhicl::ParameterSet& psphys) {

    std::string const spectrumShape = psphys.get<std::string>("spectrumShape");

    if (spectrumShape == "Bistirlich") {
      spectrum_.initialize<PionCaptureSpectrum>(elow_, ehi_, psphys.get<double>("spectrumResolution"));
    }
    else if (spectrumShape == "flat") {
      spectrum_.initialize<SimpleSpectrum>(elow_, ehi_, ehi_-elow_, SimpleSpectrum::Spectrum::Flat);
    }
    else {
      throw cet::exception("BADCONFIG")
	<< "StoppedParticlePionGun: unknown spectrum shape "<< spectrumShape << "\n";
    }
    //    return res;
  }

  //================================================================
  double StoppedParticleIntConvGun::generateEnergy() {
    return elow_ + (ehi_ - elow_)*randSpectrum_->fire();
  }

  //================================================================
  void StoppedParticleIntConvGun::produce(art::Event& event)
  {
    const auto& stop = stops_.fire();

    CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    double energy = generateEnergy();

    CLHEP::HepLorentzVector mome, momp;
    pionCaptureSpectrum_.getElecPosiVectors(energy,mome,momp); 

					// Add particles to list

    auto output = std::make_unique<GenParticleCollection>();
    output->emplace_back(PDGCode::e_minus, GenId::internalRPC,pos,mome,stop.t);
    output->emplace_back(PDGCode::e_plus , GenId::internalRPC,pos,momp,stop.t);
    event.put(move(output));
					// calculate survival probability

    double weight = exp(-stop.tauNormalized);
    event.put(std::make_unique<EventWeight>(weight));

    if (doHistograms_) {
      _hmomentum->Fill(energy);
      _hElecMom ->Fill(mome.vect().mag());
      _hPosiMom ->Fill(momp.vect().mag());

      double mee = (mome+momp).m();
      _hMee->Fill(mee);
      _hMeeVsE->Fill(energy,mee);
      _hMeeOverE->Fill(mee/energy);
      
      CLHEP::Hep3Vector p = mome.vect()+momp.vect();
      double y = (mome.e()-momp.e())/p.mag();
      
      _hy->Fill(y);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticleIntConvGun);
