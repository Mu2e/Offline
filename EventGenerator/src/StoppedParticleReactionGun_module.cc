// Andrei Gaponenko, 2013

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

#include "TH1.h"

namespace mu2e {

  //================================================================
  class StoppedParticleReactionGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    PDGCode::type       pdgId_;
    double              mass_;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };
    SpectrumVar       spectrumVariable_;

    double            elow_; // BinnedSpectrum does not store emin and emax reliably
    double            ehi_;
    BinnedSpectrum    spectrum_;
    GenId             genId_;
    int               verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    RandomUnitSphere   randomUnitSphere_;

    RootTreeSampler<IO::StoppedParticleF> stops_;
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
    double                generateEnergy();
    static BinnedSpectrum parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                             PDGCode::type              pdgId,
                                             double                     *elow,
                                             double                     *ehi);
    
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
    , ehi_ ()
    , spectrum_(parseSpectrumShape(psphys_, pdgId_, &elow_, &ehi_))
    , genId_(GenId::findByName(psphys_.get<std::string>("genId", "StoppedParticleReactionGun")))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomUnitSphere_(eng_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , doHistograms_       (pset.get<bool>("doHistograms",true ) )
  {
    produces<mu2e::GenParticleCollection>();

    if(genId_ == GenId::enum_type::unknown) {
      throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGun: unknown genId "
                                       <<psphys_.get<std::string>("genId", "StoppedParticleReactionGun")
                                       <<"\n";
    }

    if (verbosityLevel_ > 0) {
      std::cout<<"StoppedParticleReactionGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;

      std::cout<<"StoppedParticleReactionGun: using GenId = " << genId_ << std::endl;

      std::cout<<"StoppedParticleReactionGun: producing particle "<< pdgId_ << ", mass = "<< mass_ << std::endl;

      std::cout <<"StoppedParticleReactionGun: spectrum shape = "
		<<psphys_.get<std::string>("spectrumShape") << std::endl;
      if (psphys_.get<std::string>("spectrumShape")  == "tabulated")
	std::cout << " Spectrum file = "
		  << psphys_.get<std::string>("spectrumFileName")
		  << std::endl;
    }
    if (verbosityLevel_ > 1){
      std::cout <<"StoppedParticleReactionGun: spectrum: " << std::endl;
      spectrum_.print();
    }

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      //      art::TFileDirectory tfdir = tfs->mkdir( "StoppedParticleReactionGun");
      _hEnergy = tfs->make<TH1F>("hEnergy", "Energy"      , 2400,   0.0,  120);
      _hGenId  = tfs->make<TH1F>("hGenId" , "Generator ID",  100,   0.0,  100);
      _hPdgId  = tfs->make<TH1F>("hPdgId" , "PDG ID"      ,  500,  -250, 250);
      _hTime   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,   0.0, 2000.);
      _hZ      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400, 6400);

    }
  }


  //================================================================
  StoppedParticleReactionGun::SpectrumVar StoppedParticleReactionGun::parseSpectrumVar(const std::string& name) {
    if (name == "totalEnergy"  )  return TOTAL_ENERGY;
    if (name == "kineticEnergy")  return KINETIC_ENERY;
    if (name == "momentum"     )  return MOMENTUM;
    throw cet::exception("BADCONFIG")<<"StoppedParticleReactionGun: unknown spectrum variable "<<name<<"\n";
  }

  //================================================================
  BinnedSpectrum StoppedParticleReactionGun::parseSpectrumShape(const fhicl::ParameterSet& psphys,
								PDGCode::type              pdgId,
								double                     *elow,
								double                     *ehi)
  {
    BinnedSpectrum res;

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    if (spectrumShape == "Czarnecki") {
      const double mass(GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value());
      *elow = psphys.get<double>("elow", mass);
      *ehi  = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy();
      res.initialize< CzarneckiSpectrum >(*elow, *ehi, psphys.get<double>("spectrumResolution"));
    }
    else if (spectrumShape == "flat") {
      *elow = psphys.get<double>("elow");
      *ehi  = psphys.get<double>("ehi");
      res.initialize<SimpleSpectrum>(*elow, *ehi, *ehi-*elow, SimpleSpectrum::Spectrum::Flat );
    }
    else if (spectrumShape == "CeEndpoint") {
      // A simple kludge: ignore the random distribution by setting elow=ehi=eConversion
      //      *elow = *ehi = GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy();
      // default to conversion e- energy on Al
      *elow = *ehi = psphys.get<double>("ehi", GlobalConstantsHandle<PhysicsParams>()->getEndpointEnergy());
      res.initialize<SimpleSpectrum>( 0., 1., 1., SimpleSpectrum::Spectrum::Flat );
    }
    else if (spectrumShape == "ceLeadingLog") {
//-----------------------------------------------------------------------------
// ehi determines the conversion electron energy (tabulated in the .FCL file)
      //   *elow = psphys.get<double>("elow");
      *elow = 0;
      *ehi  = psphys.get<double>("ehi" );
					// for radiatively corrected spectrum, elow and ehi are derivatives 
      double bin   = psphys.get<double>("spectrumResolution");
      // int    ratio = *ehi/bin;
      // *ehi         = (ratio+1.)*bin;
      
      res.initialize<ConversionSpectrum>(*elow,*ehi,bin,*ehi,bin);
    }
    else if (spectrumShape == "ejectedProtons") {
      *elow = 0.;
      *ehi = 105.; // cut off at muon mass
      double bin = (*ehi - *elow)/psphys.get<unsigned>("nbins");
      res.initialize<EjectedProtonSpectrum>(*elow, *ehi, bin);
    }
    else if (spectrumShape == "tabulated") {
//-----------------------------------------------------------------------------
// P.Murat: assume that tabulated are the bin centers
// set 'BinCenter' to false if it is the left edges
//-----------------------------------------------------------------------------
      res.initialize(loadTable<2>( ConfigFileLookupPolicy()( psphys.get<std::string>("spectrumFileName"))) );
      *elow = res.getAbscissa(0);
      *ehi  = res.getAbscissa(res.getNbins()-1) + res.getBinWidth();
      if (psphys.get<bool>("BinCenter", false)) {
        *elow -= res.getBinWidth()/2;
        *ehi  -= res.getBinWidth()/2;
      }
      if(*elow < 0.0) throw cet::exception("BADCONFIG")
        << "StoppedParticleReactionGun: negative energy endpoint "<< *elow <<"\n";
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
                         genId_,
                         pos,
                         fourmom,
                         stop.t);

    event.put(std::move(output));
//-----------------------------------------------------------------------------
// if requested, fill histograms. Currently, the only one
//-----------------------------------------------------------------------------
    if (doHistograms_) {
      _hGenId->Fill(genId_.id());
      _hPdgId->Fill(pdgId_);
      _hEnergy->Fill(energy);
      _hTime->Fill(stop.t);
      _hZ->Fill(pos.z());
    }
  }

//-----------------------------------------------------------------------------
// generate (pseudo-)random particle energy 
// the spectrum itself doesn't know whether is stored momentum, kinetic or full 
// energy
//-----------------------------------------------------------------------------
  double StoppedParticleReactionGun::generateEnergy() {
    double res = elow_ + (ehi_ - elow_)*randSpectrum_.fire();

    if (res < 0.0) {
      throw cet::exception("BADE")<<"StoppedParticleReactionGun: negative energy "<< res <<"\n";
    }

    switch(spectrumVariable_) {
    case TOTAL_ENERGY  : break;
    case KINETIC_ENERY : res += mass_; break;
    case MOMENTUM      : res = sqrt(res*res+mass_*mass_); break;
    }
    return res;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedParticleReactionGun);
