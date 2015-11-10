// Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "DataProducts/inc/PDGCode.hh"
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
  class StoppedParticleReactionGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    PDGCode::type pdgId_;
    double mass_;

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

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    RandomUnitSphere randomUnitSphere_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    double generateEnergy();

  public:
    explicit StoppedParticleReactionGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedParticleReactionGun::StoppedParticleReactionGun(const fhicl::ParameterSet& pset)
    : psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , pdgId_(PDGCode::type(psphys_.get<int>("pdgId")))
    , mass_(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId_).ref().mass().value())
    , spectrumVariable_(parseSpectrumVar(psphys_.get<std::string>("spectrumVariable")))
    , elow_()
    , ehi_()
    , spectrum_(parseSpectrumShape(psphys_, pdgId_, &elow_, &ehi_))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomUnitSphere_(eng_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))
  {
    produces<mu2e::GenParticleCollection>();

    if(verbosityLevel_ > 0) {
      std::cout<<"StoppedParticleReactionGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StoppedParticleReactionGun: producing particle "
               <<pdgId_
               <<", mass = "<<mass_
               <<std::endl;

      std::cout <<"StoppedParticleReactionGun: spectrum shape = "
	  <<psphys_.get<std::string>("spectrumShape") << std::endl;
      if (psphys_.get<std::string>("spectrumShape")  == "tabulated")
	  std::cout << " Spectrum file = "
	  << psphys_.get<std::string>("spectrumFileName")
	  << std::endl;
    }
    if(verbosityLevel_ > 1){
      std::cout <<"StoppedParticleReactionGun: spectrum: " << std::endl;
      spectrum_.print();
    }
  }

  //================================================================
  StoppedParticleReactionGun::SpectrumVar
  StoppedParticleReactionGun::parseSpectrumVar(const std::string& name) {
    if(name == "totalEnergy")  return TOTAL_ENERGY;
    if(name == "kineticEnergy")  return KINETIC_ENERY;
    throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGun: unknown spectrum variable "<<name<<"\n";
  }

  //================================================================
  BinnedSpectrum
  StoppedParticleReactionGun::parseSpectrumShape(const fhicl::ParameterSet& psphys,
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
        << "StoppedParticleReactionGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }

    return res;
  }

  //================================================================
  void StoppedParticleReactionGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

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

    event.put(std::move(output));
  }

  //================================================================
  double StoppedParticleReactionGun::generateEnergy() {
    double res = elow_ + (ehi_ - elow_)*randSpectrum_.fire();
    switch(spectrumVariable_) {
    case TOTAL_ENERGY: break;
    case KINETIC_ENERY: res += mass_; break;
    }
    return res;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticleReactionGun);
