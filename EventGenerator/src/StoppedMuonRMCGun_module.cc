// Robert Bernstein, April 2018

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
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/Random2Dpair.hh"
#include "Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
namespace mu2e {

  //================================================================
  class StoppedMuonRMCGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    double elow_; // BinnedSpectrum does not store emin and emax reliably
    double ehi_;

    double rhoInternal_;

    BinnedSpectrum spectrum_;
    static BinnedSpectrum parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                             double *elow,
                                             double *ehi);

    int verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    RandomUnitSphere randomUnitSphere_;

    CLHEP::RandFlat randFlat_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    bool kMaxUserSet_;
    double kMaxUser_;
    bool doHistograms_;
 
    static double fractionSpectrum ;
    static double omcNormalization;

    double generateEnergy();

    TH1F* _hmomentum;
    TH1F* _hEnergyElectron;
    TH1F* _hEnergyPositron;
    TH1F* _hWeight;

  public:
    explicit StoppedMuonRMCGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  double StoppedMuonRMCGun::fractionSpectrum{0.};
  double StoppedMuonRMCGun::omcNormalization{0.};
  //  double StoppedMuonRMCGun::rhoInternal{0.0069};

  //================================================================
  StoppedMuonRMCGun::StoppedMuonRMCGun(const fhicl::ParameterSet& pset)
    : psphys_(pset.get<fhicl::ParameterSet>("physics"))
    , elow_()
    , ehi_()
    , rhoInternal_(psphys_.get<double>("rhoInternal"))
    , spectrum_(parseSpectrumShape(psphys_, &elow_, &ehi_))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_(eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randomUnitSphere_(eng_)
    , randFlat_(eng_)
    , stops_(eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , doHistograms_( pset.get<bool>("doHistograms",true ) )
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::EventWeight>();

    if(verbosityLevel_ > 0) {
      std::cout<<"StoppedMuonRMCGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;
      std::cout<<"StoppedMuonRMCGun: producing photon " << std::endl;
    }



    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "StoppedMuonRMCGun" );

      _hmomentum     = tfdir.make<TH1F>( "hmomentum", "Produced photon momentum, RMC", 70,  0.,  140.  );
      _hEnergyElectron     = tfdir.make<TH1F>( "hEnergyElectron", "Produced electron energy, RMC Internal", 70,  0.,  140.  );
      _hEnergyPositron     = tfdir.make<TH1F>( "hEnergyPositron", "Produced electron energy, RMC Internal", 70,  0.,  140.  );
      _hWeight             = tfdir.make<TH1F>( "hWeight",         "Event Weight ", 100,0.,1.);
    }

  }

  //================================================================
  BinnedSpectrum 
  StoppedMuonRMCGun::parseSpectrumShape(const fhicl::ParameterSet& psphys,
                                                 double *elow,
                                                 double *ehi)
  {
    BinnedSpectrum res;

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    const int physicsVerbosityLevel_(psphys.get<int>("physicsVerbosityLevel"));
    if (spectrumShape == "ClosureApprox") {
      bool blind = psphys.get<bool>("blind");
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      bool kMaxUserSet = psphys.get<bool>("kMaxUserSet");
      double kMaxUser = psphys.get<double>("kMaxUser");
      double rmcFrac =  psphys.get<double>("rmcFrac");
 

    // in this case just stop
      if (*elow >= *ehi){
	// this is wrong
	throw cet::exception("RANGE") << "energy range in Muon Capture Spectrum is wrong " << *elow << " " << *ehi << std::endl;
      }

      if (physicsVerbosityLevel_ > 0 && !blind){
	std::cout << "kMaxUserSet and kMaxUser = " << kMaxUserSet << " " << kMaxUser << std::endl;
      }
      //
      // this will allow me to select a region of the RMC spectrum and weight for the fraction of the spectrum
      const double muonMassFit{105.658};
      const double bindingEnergyFit{0.464};
      const double recoilEnergyFit{0.220};
      const double deltaMassFit{3.121};
      const double kMaxMax{muonMassFit - bindingEnergyFit - recoilEnergyFit - deltaMassFit};
 
      double kMax;
      if (kMaxUserSet){
	kMax = kMaxUser;
      } 
      else {
	kMax = kMaxMax;
      } 
 
      if (  *elow > kMax ) {
	//
	// if I told you what kMax was you could unblind kMax.  Therefore I will set it to something very low and tell you.
	std::cout << " StoppedMuonGun elow is too high " << *elow << " resetting to 0 MeV" << std::endl;

      }

      res.initialize<MuonCaptureSpectrum>( *elow, *ehi, psphys.get<double>("spectrumResolution"),
					   kMaxUserSet, kMaxUser, kMaxMax);
 
 
 

      double lowestEnergy{*elow};
      double upperEnergy{*ehi};

      if (*ehi > kMax) {
	upperEnergy = kMax;
      }

 
      double xLower = lowestEnergy/kMax;
      double xUpper = upperEnergy/kMax;
    
      //
      // integral of closure appoximation is 1/20 over [0,1].  this gives me the fraction of the spectrum we use, 
      // should weight overall rate by this

      fractionSpectrum = (20.) *  ( pow(xUpper,2)/2. - (4./3.)*pow(xUpper,3) + (7./4.)*pow(xUpper,4) - (6./5.)*pow(xUpper,5) 
			   + (1./3.)*pow(xUpper,6) )
	-  
	( pow(xLower,2)/2. - (4./3.)*pow(xLower,3) + (7./4.)*pow(xLower,4) - (6./5.)*pow(xLower,5) 
	  + (1./3.)*pow(xLower,6) );

      //
      // this is a DIFFERENT normalization.  Docdb 4378 and Armstrong et al tell us the rate about 57 MeV normalized to all
      // ordinary muon captures is 1.43 +-0.12 x 10^{-5}.  See the mathematica notebook in doc-db 16979.  Made configurable.

      omcNormalization = (rmcFrac)/(1/20. - 11432149083/pow(kMax,6) + (3610152342/5.)/pow(kMax,5) - (73892007/4.)/pow(kMax,4) + 246924/pow(kMax,3) - (3249/2.)/pow(kMax,2));
 
      if (physicsVerbosityLevel_ > 0){
	std::cout << "fraction of spectrum = " << fractionSpectrum << std::endl;
	std::cout << "rmc fraction         = " << rmcFrac << std::endl;
	std::cout << "omc normalization    = " << omcNormalization << std::endl;
      }

    }
    else if (spectrumShape == "flat") {
      *elow = psphys.get<double>("elow");
      *ehi = psphys.get<double>("ehi");
      res.initialize<SimpleSpectrum>(*elow, *ehi, *ehi-*elow, SimpleSpectrum::Spectrum::Flat );
    }
    else {
      throw cet::exception("BADCONFIG")
        << "StoppedParticleMuonGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }

    return res;
  }

  //================================================================
  void StoppedMuonRMCGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    const double energy = generateEnergy();

    double weight{0.};
    // two things can now happen with this photon.  It can proceed and possibly convert, or it can internally convert.
    // the probability of internal conversion is given by rho = 0.0069 (assume same as for RPC).  Throw a flat 

    if (randFlat_.fire() > rhoInternal_) {
      output->emplace_back( PDGCode::gamma, 
			    GenId::radiativeMuonCapture, 
			    pos,
			    CLHEP::HepLorentzVector( randomUnitSphere_.fire(energy), energy), 
			    stop.t );

      event.put(std::move(output));

      // for future normalization
      const double weightExternal = fractionSpectrum * omcNormalization;
      std::unique_ptr<EventWeight> pw(new EventWeight(weightExternal));
      event.put(std::move(pw));
      weight = weightExternal;
      if ( doHistograms_ ) {

	_hmomentum->Fill(energy);
      }
    } else {    // internal conversions

      //Need mass of electron
      static const double massE = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus).ref().mass().value();

      // Uses energy above as photon energy, assuming distribution created by 
      static Random2Dpair< MuonCaptureSpectrum > random2dPair( eng_, 2*massE, energy, -1., 1. );
 
      const auto xyPair          = random2dPair.fire( energy );
      const auto elecPosiVectors = MuonCaptureSpectrum::getElecPosiVectors( energy, xyPair.first, xyPair.second ); 
      CLHEP::HepLorentzVector fakeElectron( 105.*TMath::Sin(CLHEP::pi*60./180.),0.,105.*TMath::Cos(CLHEP::pi*60./180.),sqrt(105*105+ massE*massE));
      CLHEP::HepLorentzVector fakePositron(-105.*TMath::Sin(CLHEP::pi*60./180.),0.,105.*TMath::Cos(CLHEP::pi*60./180.),sqrt(105*105+ massE*massE));
      output->emplace_back( PDGCode::e_minus, 
			    GenId::radiativeMuonCaptureInternal, 
			    pos,
			    //			    elecPosiVectors.first, 
			    fakeElectron, 
			    800. );
      //			    stop.t );
      output->emplace_back( PDGCode::e_plus, 
			    GenId::radiativeMuonCaptureInternal, 
			    pos,
			    //			    elecPosiVectors.second, 
			    fakePositron, 
			    800.);
      //			    stop.t );

      event.put(std::move(output));

      // for future normalization
      const double weightInternal = omcNormalization*fractionSpectrum*rhoInternal_;
      std::unique_ptr<EventWeight> pw(new EventWeight(weightInternal));
      weight = weightInternal;
      event.put(std::move(pw));
      if (verbosityLevel_ > 0) {
	std::cout << "original photon energy = " << energy << " and electron mass = " << massE <<  std::endl;
	std::cout << "RMC electron/positron energies = " << elecPosiVectors.first.e() << " " << elecPosiVectors.second.e() << std::endl;
	std::cout << "and the full 4-vector: " << elecPosiVectors.first << " " << elecPosiVectors.second << std::endl;
	std::cout << "stop time = " << stop.t << std::endl;
	std::cout << " event weight = " << fractionSpectrum << " " << rhoInternal_ << " " << weightInternal << std::endl;
      }

      if ( doHistograms_ ) {
	_hWeight->Fill(weight);
	_hmomentum->Fill(energy);
	_hEnergyElectron->Fill(elecPosiVectors.first.e());
	_hEnergyPositron->Fill(elecPosiVectors.second.e());

      }
    }
  }

  //================================================================
  double StoppedMuonRMCGun::generateEnergy() {

    return elow_ + (ehi_ - elow_)*randSpectrum_.fire();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedMuonRMCGun);
