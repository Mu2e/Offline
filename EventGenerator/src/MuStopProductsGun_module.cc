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
#include "CLHEP/Random/RandPoissonQ.h"
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

      fhicl::Sequence< fhicl::Table<GenPhysConfig> > stopProducts{Name("stopProducts"), Comment("Products coincident with stopped muon (i.e. generated every event)")};
      fhicl::Sequence< fhicl::Table<GenPhysConfig> > captureProducts{Name("captureProducts"), Comment("Products coincident with captured muon)")};
      fhicl::Sequence< fhicl::Table<GenPhysConfig> > decayProducts{Name("decayProducts"), Comment("Products coincident with decayd muon)")};
      fhicl::Atom<int> verbosityLevel{Name("verbosityLevel"), Comment("Verbosity Level (default = 0)"), 0};
      fhicl::Table<RTS::Config> stops{Name("stops"), Comment("Stops ntuple config")};
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), Comment("True/false to produce histograms"), true};
      fhicl::Atom<bool> doTree{Name("doTree"), Comment("True/false to produce tree"), false};
    };
    typedef art::EDProducer::Table<Config> Parameters;

  private:
    Config conf_;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };

    struct GenPhysStruct {
      GenPhysStruct(GenPhysConfig conf, art::RandomNumberGenerator::base_engine_t& eng, bool is_poisson_rate = false, double rate = 1, double prob = 1) : 
	pdgId(PDGCode::type(conf.pdgId())),
	mass(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value()),
	spectrumVariable(parseSpectrumVar(conf.spectrumVariable())),
	spectrum(BinnedSpectrum(conf)),
	genId(GenId::findByName(conf.genId())),
	randSpectrum(eng, spectrum.getPDF(), spectrum.getNbins()),
	isPoissonRate(is_poisson_rate),
	emissionRate(rate),
	probability(prob),
	randomFlat(eng),
	randomPoissonQ(eng,rate)
      {    }

      PDGCode::type       pdgId;
      double              mass;
      SpectrumVar       spectrumVariable;
      BinnedSpectrum    spectrum;
      GenId             genId;
      CLHEP::RandGeneral randSpectrum;
      bool isPoissonRate;
      double emissionRate;
      double probability;
      CLHEP::RandFlat randomFlat;
      CLHEP::RandPoissonQ randomPoissonQ;

      int getNGenerated() {
	double rand = randomFlat.fire();
	if (rand < probability) {
	  if (isPoissonRate) {
	    return randomPoissonQ.fire();
	  }
	  else {
	    return emissionRate;
	  }
	}
	else {
	  return 0;
	}
      };
    };


    int               verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    RandomUnitSphere   randomUnitSphere_;
    CLHEP::RandFlat randomFlat_;

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
    TH1F*   _hProcesses;

    bool    doTree_;
    TTree*  _genTree;
    Float_t _brEnergy;
    Int_t   _brGenId;
    Int_t   _brPdgId;
    Float_t _brTime;
    Float_t _brZ;
  private:
    static SpectrumVar    parseSpectrumVar(const std::string& name);
    double                generateEnergy(GenPhysStruct& genPhys);
    void                  generateGenParticles(GenPhysStruct& genPhys, std::unique_ptr<GenParticleCollection>& output);
    double _decayFraction;
    double _captureFraction;
    std::vector<GenPhysStruct> _allStopProducts;
    std::vector<GenPhysStruct> _allCaptureProducts;
    std::vector<GenPhysStruct> _allDecayProducts;

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
    , randomFlat_(eng_)
    , stops_(eng_, conf_.stops())
    , doHistograms_(conf_.doHistograms())
    , doTree_(conf_.doTree())
    , _decayFraction(GlobalConstantsHandle<PhysicsParams>()->getDecayFraction())
    , _captureFraction(1 - _decayFraction)
  {
    produces<mu2e::GenParticleCollection>();

    if (verbosityLevel_ > 0) {
      std::cout<<"MuStopProductsGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;
      std::cout << "MuStopProductsGun: decayFraction = " << _decayFraction << std::endl;
      std::cout << "MuStopProductsGun: captureFraction = " << _captureFraction << std::endl;
    }
    for (const auto& i_genPhysConfig : conf_.stopProducts()) {
      double rate = 1;
      bool is_poisson_rate = false;
      double prob = 1;
      if (i_genPhysConfig.pdgId() == 22) {
	prob = 0.79;//rate = GlobalConstantsHandle<PhysicsParams>()->getCaptureProtonRate();
      }

      GenPhysStruct i_genPhys(i_genPhysConfig, eng_, is_poisson_rate, rate, prob);
      _allStopProducts.push_back(i_genPhys);
      if (verbosityLevel_ > 0) {
	std::cout << "MuStopProductsGun: Adding Stop Process: " << i_genPhys.pdgId << std::endl;
      }
    }
    for (const auto& i_genPhysConfig : conf_.captureProducts()) {
      double rate = 0;
      bool is_poisson_rate = true;
      double prob = 1;
      if (i_genPhysConfig.pdgId() == 2212) {
	rate = GlobalConstantsHandle<PhysicsParams>()->getCaptureProtonRate();
      }
      else if (i_genPhysConfig.pdgId() == 1000010020) {
	rate = GlobalConstantsHandle<PhysicsParams>()->getCaptureDeuteronRate();
      }
      else if (i_genPhysConfig.pdgId() == 2112) {
	rate = GlobalConstantsHandle<PhysicsParams>()->getCaptureNeutronRate();
      }
      else {
	throw cet::exception("MUSTOPPRODUCTSGUN") << "Capture product with PdgId " << i_genPhysConfig.pdgId() << " is not implemented" << std::endl;
      }

      GenPhysStruct i_genPhys(i_genPhysConfig, eng_, is_poisson_rate, rate, prob);
      _allCaptureProducts.push_back(i_genPhys);
      if (verbosityLevel_ > 0) {
	std::cout << "MuStopProductsGun: Adding Capture Process: " << i_genPhys.pdgId << std::endl;
      }
    }
    for (const auto& i_genPhysConfig : conf_.decayProducts()) {
      double rate = 0;
      bool is_poisson_rate = false;
      double prob = 1;
      if (i_genPhysConfig.pdgId() == 11) {
	rate = 1;
      }
      else {
	throw cet::exception("MUSTOPPRODUCTSGUN") << "Emission rate/Probability for decay product with PdgId " << i_genPhysConfig.pdgId() << " is not implemented" << std::endl;
      }

      GenPhysStruct i_genPhys(i_genPhysConfig, eng_, is_poisson_rate, rate, prob);
      _allDecayProducts.push_back(i_genPhys);
      if (verbosityLevel_ > 0) {
	std::cout << "MuStopProductsGun: Adding Decay Process: " << i_genPhys.pdgId << std::endl;
      }
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      //      art::TFileDirectory tfdir = tfs->mkdir( "MuStopProductsGun");
      _hEnergy = tfs->make<TH1F>("hEnergy", "Energy"      , 2400,   0.0,  120);
      _hGenId  = tfs->make<TH1F>("hGenId" , "Generator ID",  100,   0.0,  100);
      _hPdgId  = tfs->make<TH1F>("hPdgId" , "PDG ID"      ,  500,  -250, 250);
      _hTime   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,   0.0, 2000.);
      _hZ      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400, 6400);

      _hProcesses = tfs->make<TH1F>("hProcesses", "Process Counts", 3,0,3);
      _hProcesses->GetXaxis()->SetBinLabel(1, "stops");
      _hProcesses->GetXaxis()->SetBinLabel(2, "decays");
      _hProcesses->GetXaxis()->SetBinLabel(3, "captures");
    }
    if ( doTree_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      //      art::TFileDirectory tfdir = tfs->mkdir( "MuStopProductsGun");
      _genTree = tfs->make<TTree>("genTree", "tree of generated particles");
      _genTree->Branch("energy", &_brEnergy);
      _genTree->Branch("genId", &_brGenId);
      _genTree->Branch("pdgId", &_brPdgId);
      _genTree->Branch("time", &_brTime);
      _genTree->Branch("z", &_brZ);
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

    // Always generate stop products (e.g. stop X-rays)
    if (doHistograms_) {
      _hProcesses->Fill("stops", 1);
    }
    for (auto& i_genPhys : _allStopProducts) {
      generateGenParticles(i_genPhys, output);
    }

    double rand = randomFlat_.fire();
    if (rand < _decayFraction) {
      if (doHistograms_) {
	_hProcesses->Fill("decays", 1);
      }
      for (auto& i_genPhys : _allDecayProducts) {
	generateGenParticles(i_genPhys, output);
      }
    }
    else {
      if (doHistograms_) {
	_hProcesses->Fill("captures", 1);
      }
      for (auto& i_genPhys : _allCaptureProducts) {
	generateGenParticles(i_genPhys, output);
      }
    }

    event.put(std::move(output));
  }

//-----------------------------------------------------------------------------
// generate (pseudo-)random particle energy 
// the spectrum itself doesn't know whether is stored momentum, kinetic or full 
// energy
//-----------------------------------------------------------------------------
  double MuStopProductsGun::generateEnergy(GenPhysStruct& genPhys) {
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

  // Generate the GenParticles and add them to the output
  void MuStopProductsGun::generateGenParticles(GenPhysStruct& genPhys, std::unique_ptr<GenParticleCollection>& output) {

    const auto& stop = stops_.fire();
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    int n_gen = genPhys.getNGenerated();

    for (int i_gen = 0; i_gen < n_gen; ++i_gen) {
      const double energy = generateEnergy(genPhys);
      const double p = energy * sqrt(1 - std::pow(genPhys.mass/energy,2));

      CLHEP::Hep3Vector p3 = randomUnitSphere_.fire(p);
      CLHEP::HepLorentzVector fourmom(p3, energy);
      output->emplace_back(genPhys.pdgId,
			   genPhys.genId,
			   pos,
			   fourmom,
			   stop.t);

      //-----------------------------------------------------------------------------
      // if requested, fill histograms. Currently, the only one
      //-----------------------------------------------------------------------------
      if (doHistograms_) {
	_hGenId->Fill(genPhys.genId.id());
	_hPdgId->Fill(genPhys.pdgId);
	_hEnergy->Fill(energy);
	_hTime->Fill(stop.t);
	_hZ->Fill(pos.z());
      }
      if (doTree_) {
	_brEnergy = energy;
	_brGenId = genPhys.genId.id();
	_brPdgId = genPhys.pdgId;
	_brTime = stop.t;
	_brZ = pos.z();

	_genTree->Fill();
	
	// Reset branches
	_brEnergy = 0;
	_brGenId = 0;
	_brPdgId = 0;
	_brTime = 0;
	_brZ = 0;
      }
    }
  }
  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::MuStopProductsGun);
