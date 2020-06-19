// Andy Edmonds, 2020
// based on StoppedParticleReactionGun by Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

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

#include "fhiclcpp/types/OptionalAtom.h"

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
#include "Mu2eUtilities/inc/ConversionSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "Mu2eUtilities/inc/GenPhysConfig.hh"

#include "TH1.h"

namespace mu2e {

  //================================================================
  class MuStopProductsGun : public art::EDProducer {

    typedef RootTreeSampler<IO::StoppedParticleF> RTS;

  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence< fhicl::Table<GenPhysConfig> > physics{Name("physics"), Comment("Physics parameter fhicl table")};
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity Level (default = 0)"), 0};
      fhicl::Table<RTS::Config> stops{Name("stops"), Comment("Stops ntuple config")};
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), Comment("True/false to produce histograms"), true};
    };
    typedef art::EDProducer::Table<Config> Parameters;

  private:
    Config conf_;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };

    struct GenPhys {
      GenPhys(GenPhysConfig conf, art::RandomNumberGenerator::base_engine_t& eng) : 
	pdgId(PDGCode::type(conf.pdgId())),
	mass(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value()),
	spectrumVariable(parseSpectrumVar(conf.spectrumVariable())),
	spectrum(BinnedSpectrum(conf)),
	genId(GenId::findByName(conf.genId())),
	randSpectrum(eng, spectrum.getPDF(), spectrum.getNbins())
      {}

      PDGCode::type       pdgId;
      double              mass;
      SpectrumVar       spectrumVariable;
      BinnedSpectrum    spectrum;
      GenId             genId;
      CLHEP::RandGeneral randSpectrum;
    };
    std::vector<GenPhys> allGenPhys_;

    int               verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    RandomUnitSphere   randomUnitSphere_;

    RTS stops_;
//-----------------------------------------------------------------------------
// histogramming
//-----------------------------------------------------------------------------
    bool    doHistograms_;
    TH1F*   _hEnergy;
    TH1F*   _hPdgId;
    TH1F*   _hGenId;
    TH1F*   _hTime;
    TH1F*   _hZ;
  
  private:
    static SpectrumVar    parseSpectrumVar(const std::string& name);
    double                generateEnergy(GenPhys& genPhys);
    
  public:
    explicit MuStopProductsGun(const Parameters& conf);

    virtual void produce(art::Event& event);
  };

  //================================================================
  MuStopProductsGun::MuStopProductsGun(const Parameters& conf)
    : EDProducer(conf)
    , conf_(conf())
    , verbosityLevel_(conf_.verbosityLevel())
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randomUnitSphere_(eng_)
    , stops_(eng_, conf_.stops())
    , doHistograms_(conf_.doHistograms())
  {
    produces<mu2e::GenParticleCollection>();

    for (const auto& i_genPhysConfig : conf_.physics()) {
      GenPhys i_genPhys(i_genPhysConfig, eng_);

      if(i_genPhys.genId == GenId::enum_type::unknown) {
	throw cet::exception("BADCONFIG")<<"MuStopProductsGun: unknown genId "
					 << i_genPhysConfig.genId()
					 <<"\n";
      }

      if (verbosityLevel_ > 0) {
	std::cout<<"MuStopProductsGun: using GenId = " << i_genPhys.genId << std::endl;

	std::cout<<"MuStopProductsGun: producing particle "<< i_genPhys.pdgId << ", mass = "<< i_genPhys.mass << std::endl;

	std::cout <<"MuStopProductsGun: spectrum shape = "
		  << i_genPhysConfig.spectrumShape() << std::endl;
	if (i_genPhysConfig.spectrumShape()  == "tabulated") {
	  std::string spectrumFileName;
	  if (i_genPhysConfig.spectrumFileName(spectrumFileName)) {
	    std::cout << " Spectrum file = "
		      << spectrumFileName
		      << std::endl;
	  }
	}
      }
      if (verbosityLevel_ > 1){
	std::cout <<"MuStopProductsGun: spectrum: " << std::endl;
	i_genPhys.spectrum.print();
      }

      allGenPhys_.push_back(i_genPhys);
    }

    if (verbosityLevel_ > 0) {
      std::cout<<"MuStopProductsGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      //      art::TFileDirectory tfdir = tfs->mkdir( "MuStopProductsGun");
      _hEnergy = tfs->make<TH1F>("hEnergy", "Energy"      , 2400,   0.0,  120);
      _hGenId  = tfs->make<TH1F>("hGenId" , "Generator ID",  100,   0.0,  100);
      _hPdgId  = tfs->make<TH1F>("hPdgId" , "PDG ID"      ,  500,  -250, 250);
      _hTime   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,   0.0, 2000.);
      _hZ      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400, 6400);

    }
  }


  //================================================================
  MuStopProductsGun::SpectrumVar MuStopProductsGun::parseSpectrumVar(const std::string& name) {
    if (name == "totalEnergy"  )  return TOTAL_ENERGY;
    if (name == "kineticEnergy")  return KINETIC_ENERY;
    if (name == "momentum"     )  return MOMENTUM;
    throw cet::exception("BADCONFIG")<<"MuStopProductsGun: unknown spectrum variable "<<name<<"\n";
  }


  //================================================================
  void MuStopProductsGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    for (auto& i_genPhys : allGenPhys_) {
      const double energy = generateEnergy(i_genPhys);
      const double p = energy * sqrt(1 - std::pow(i_genPhys.mass/energy,2));

      CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);
      output->emplace_back(i_genPhys.pdgId,
			   i_genPhys.genId,
			   pos,
			   fourmom,
			   stop.t);

      //-----------------------------------------------------------------------------
      // if requested, fill histograms. Currently, the only one
      //-----------------------------------------------------------------------------
      if (doHistograms_) {
	_hGenId->Fill(i_genPhys.genId.id());
	_hPdgId->Fill(i_genPhys.pdgId);
	_hEnergy->Fill(energy);
	_hTime->Fill(stop.t);
	_hZ->Fill(pos.z());
      }

    }
    event.put(std::move(output));
  }

//-----------------------------------------------------------------------------
// generate (pseudo-)random particle energy 
// the spectrum itself doesn't know whether is stored momentum, kinetic or full 
// energy
//-----------------------------------------------------------------------------
  double MuStopProductsGun::generateEnergy(GenPhys& genPhys) {
    double res = genPhys.spectrum.sample(genPhys.randSpectrum.fire());

    if (res < 0.0) {
      throw cet::exception("BADE")<<"MuStopProductsGun: negative energy "<< res <<"\n";
    }

    switch(genPhys.spectrumVariable) {
    case TOTAL_ENERGY  : break;
    case KINETIC_ENERY : res += genPhys.mass; break;
    case MOMENTUM      : res = sqrt(res*res+genPhys.mass*genPhys.mass); break;
    }
    return res;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuStopProductsGun);
