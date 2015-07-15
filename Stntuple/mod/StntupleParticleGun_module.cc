// Andrei Gaponenko, 2013
// cloned version which can handle multiple particles

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib/exception.h"

#include "TH1.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  //================================================================
  class StntupleParticleGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    double        mean_;
    PDGCode::type pdgId_;
    double        mass_;

    enum SpectrumVar { TOTAL_ENERGY, KINETIC_ENERY };
    SpectrumVar spectrumVariable_;
    static SpectrumVar parseSpectrumVar(const std::string& name);

    double elow_; // BinnedSpectrum does not store emin and emax reliably
    double ehi_;
    BinnedSpectrum spectrum_;
    static BinnedSpectrum parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                             PDGCode::type pdgId,
                                             double *elow,
                                             double *ehi);

    int verbosityLevel_;
    int diagLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral    randSpectrum_;
    RandomUnitSphere      randomUnitSphere_;

    CLHEP::RandPoissonQ   randPoissonQ_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    TH1F   *hEnergy_;

    double generateEnergy();

  public:
    explicit StntupleParticleGun(const fhicl::ParameterSet& pset);

    virtual void produce(art::Event& event);
  };

  //================================================================
  StntupleParticleGun::StntupleParticleGun(const fhicl::ParameterSet& pset)
    : psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , mean_  (pset.get<double>("mean",1.))
    , pdgId_(PDGCode::type(psphys_.get<int>("pdgId")))
    , mass_(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , spectrumVariable_(parseSpectrumVar(psphys_.get<std::string>("spectrumVariable")))
    , elow_()
    , ehi_()
    , spectrum_(parseSpectrumShape(psphys_, pdgId_, &elow_, &ehi_))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , diagLevel_(pset.get<int>("diagLevel", 0))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomUnitSphere_(eng_)
    , randPoissonQ_(eng_, std::abs(mean_))
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosityLevel_ > 0) {
      std::cout<<"StntupleParticleGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StntupleParticleGun: producing particle "
               <<pdgId_
               <<", mass = "<<mass_
               <<std::endl;

      std::cout <<"StntupleParticleGun: spectrum shape = "
	  <<psphys_.get<std::string>("spectrumShape") << std::endl;
      if (psphys_.get<std::string>("spectrumShape")  == "tabulated")
	  std::cout << " Spectrum file = "
	  << psphys_.get<std::string>("spectrumFileName")
	  << std::endl;
    }
    if(verbosityLevel_ > 1){
      std::cout <<"StntupleParticleGun: spectrum: " << std::endl;
      spectrum_.print();
    }

    art::ServiceHandle<art::TFileService> tfs;

    if (diagLevel_ > 0) {
        hEnergy_ = tfs->make<TH1F>("energy","Particle Energy",100,0,200);
    }
  }

  //================================================================
  StntupleParticleGun::SpectrumVar
  StntupleParticleGun::parseSpectrumVar(const std::string& name) {
    if(name == "totalEnergy")  return TOTAL_ENERGY;
    if(name == "kineticEnergy")  return KINETIC_ENERY;
    throw cet::exception("BADCONFIG")<<"StntupleParticleGun: unknown spectrum variable "<<name<<"\n";
  }

  //================================================================
  BinnedSpectrum
  StntupleParticleGun::parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                                 PDGCode::type pdgId,
                                                 double *elow,
                                                 double *ehi)
  {
    BinnedSpectrum res;

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    if (spectrumShape == "Czarnecki") {
      const double mass(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value());
      *elow = psphys.get<double>("elow", mass);
      *ehi = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy();
      res.initialize< CzarneckiSpectrum >(*elow, *ehi, psphys.get<double>("spectrumResolution"));
    }
    else if (spectrumShape == "flat") {
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      res.initialize<SimpleSpectrum>(*elow, *ehi, *ehi-*elow, SimpleSpectrum::Spectrum::Flat );
    }
    else if (spectrumShape == "conversion") {
      // A simple kludge: ignore the random distribution by setting elow=ehi=eConversion
      res.initialize<SimpleSpectrum>(0., 1., 1., SimpleSpectrum::Spectrum::Flat );
      *elow = *ehi = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy();
    }
    else if (spectrumShape == "ejectedProtons") {
      *elow = 0.;
      *ehi = 105.; // cut off at muon mass
      double spectrumRes = (*ehi - *elow)/psphys.get<unsigned>("nbins");
      res.initialize<EjectedProtonSpectrum>(*elow, *ehi, spectrumRes);
    }
    else if (spectrumShape == "tabulated") {
      res.initialize(loadTable<2>( ConfigFileLookupPolicy()( psphys.get<std::string>("spectrumFileName"))) );
      *elow = res.getAbscissa(0);
      *ehi = res.getAbscissa(res.getNbins()-1);
    }
    else {
      throw cet::exception("BADCONFIG")
        << "StntupleParticleGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }

    return res;
  }

  //================================================================
  void StntupleParticleGun::produce(art::Event& event) {

    // Choose the number of electrons to generate this event.

    long n = (mean_ < 0 ? static_cast<long>(-mean_): randPoissonQ_.fire());

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    for (int i=0; i<n; ++i) {
      const auto& stop = stops_.fire();

      const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

      const double energy = generateEnergy();
      const double p = energy * sqrt(1 - std::pow(mass_/energy,2));

      CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);

      output->emplace_back(pdgId_,
			   GenId::StoppedParticleReactionGun,
			   pos,
			   fourmom,
			   stop.t);

      if (diagLevel_ > 0) {
        hEnergy_->Fill(energy);
      }
    }

    event.put(std::move(output));
  }

  //================================================================
  double StntupleParticleGun::generateEnergy() {
    double res = elow_ + (ehi_ - elow_)*randSpectrum_.fire();
    switch(spectrumVariable_) {
    case TOTAL_ENERGY: break;
    case KINETIC_ENERY: res += mass_; break;
    }
    return res;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StntupleParticleGun);
