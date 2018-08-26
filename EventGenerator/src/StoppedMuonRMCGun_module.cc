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
// #include "Mu2eUtilities/inc/Random2Dpair.hh"
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
#include "TH2F.h"
#include "TMath.h"
namespace mu2e {

  //================================================================
  class StoppedMuonRMCGun : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    double rhoInternal_;
    double elow_; // BinnedSpectrum does not store emin and emax reliably
    double ehi_;


    BinnedSpectrum spectrum_;

    int verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGeneral randSpectrum_;
    RandomUnitSphere   randUnitSphere_;
    CLHEP::RandFlat    randFlat_;

    MuonCaptureSpectrum muonCaptureSpectrum_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    bool   kMaxUserSet_;
    double kMaxUser_;
    bool   doHistograms_;

    double me_;                        // electron mass
    double mmu_;
 
    double fractionSpectrum_;
    double omcNormalization_;

    TH1F* _hmomentum;
    TH1F* _hEnergyElectron;
    TH1F* _hEnergyPositron;
    TH1F* _hWeight;
    TH1F* _htZero;
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hy;				// splitting function
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    double generateEnergy();
    void   parseSpectrumShape(const fhicl::ParameterSet& psphys);

  public:
    explicit StoppedMuonRMCGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };


//================================================================
  StoppedMuonRMCGun::StoppedMuonRMCGun(const fhicl::ParameterSet& pset)
    : psphys_             (pset.get<fhicl::ParameterSet>("physics"))
    , rhoInternal_        (psphys_.get<double>("rhoInternal"))
    , elow_               (psphys_.get<double>("elow"))
    , ehi_                (psphys_.get<double>("ehi"))
    , verbosityLevel_     (pset.get<int>("verbosityLevel", 0))
    , eng_                (createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randSpectrum_       (eng_, spectrum_.getPDF(), spectrum_.getNbins())
    , randUnitSphere_     (eng_)
    , randFlat_           (eng_)
    , muonCaptureSpectrum_(&randFlat_,&randUnitSphere_)
    , stops_              (eng_, pset.get<fhicl::ParameterSet>("muonStops"))
    , doHistograms_       (pset.get<bool>("doHistograms",true ) )
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::EventWeight>();
    
    fractionSpectrum_ = 0.;
    omcNormalization_ = 0.;

    if(verbosityLevel_ > 0) {
      std::cout<<"StoppedMuonRMCGun: using = "
               <<stops_.numRecords()
               <<" stopped particles"
               <<std::endl;
      std::cout<<"StoppedMuonRMCGun: producing photon " << std::endl;
    }

    me_  = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus ).ref().mass().value();
    mmu_ = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::mu_minus).ref().mass().value();

    // initialize binned spectrum - this needs to be done right
    parseSpectrumShape(psphys_);

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "StoppedMuonRMCGun" );

      _hmomentum     = tfdir.make<TH1F>( "hmomentum", "Produced photon momentum, RMC", 70,  0.,  140.  );
      _hEnergyElectron     = tfdir.make<TH1F>( "hEnergyElectron", "Produced electron energy, RMC Internal", 70,  0.,  140.  );
      _hEnergyPositron     = tfdir.make<TH1F>( "hEnergyPositron", "Produced electron energy, RMC Internal", 70,  0.,  140.  );
      _htZero              = tfdir.make<TH1F>( "htZero"         , "Stopped Muon time", 100,0.,2000.);
      _hWeight             = tfdir.make<TH1F>( "hWeight"        , "Event Weight ", 100,0.,1.);
      _hMee                = tfdir.make<TH1F>( "hMee"           , "M(e+e-) "     , 200,0.,200.);
      _hMeeVsE             = tfdir.make<TH2F>( "hMeeVsE"        , "M(e+e-) "     , 200,0.,200.,200,0,200);
      _hy                  = tfdir.make<TH1F>( "hy"             , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
    }
  }

  //================================================================
  void StoppedMuonRMCGun::parseSpectrumShape(const fhicl::ParameterSet& psphys) {

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    const int physicsVerbosityLevel_(psphys.get<int>("physicsVerbosityLevel"));

    if (spectrumShape == "ClosureApprox") {
					// in this case just stop, this is wrong
      if (elow_ >= ehi_){
	throw cet::exception("RANGE") << "energy range in Muon Capture Spectrum is wrong " << elow_ << " " << ehi_ << std::endl;
      }

      bool   blind       = psphys.get<bool>  ("blind");
      bool   kMaxUserSet = psphys.get<bool>  ("kMaxUserSet");

      double kMaxUser(0);
      if (kMaxUserSet) kMaxUser = psphys.get<double>("kMaxUser");

      double rmcFrac     = psphys.get<double>("rmcFrac");
      double bin         = psphys.get<double>("spectrumResolution");
 
      if (physicsVerbosityLevel_ > 0 && !blind) {
	std::cout << "kMaxUserSet and kMaxUser = " << kMaxUserSet << " " << kMaxUser << std::endl;
      }

      //
      // this will allow me to select a region of the RMC spectrum and weight for the fraction of the spectrum

      const double bindingEnergyFit{0.464};
      const double recoilEnergyFit {0.220};
      const double deltaMassFit    {3.121};
      const double kMaxMax         {mmu_ - bindingEnergyFit - recoilEnergyFit - deltaMassFit};
 
      double kMax;
      if (kMaxUserSet) kMax = kMaxUser;
      else             kMax = kMaxMax;
 
      if ( elow_ > kMax ) {
	//
	// if I told you what kMax was you could unblind kMax.  Therefore I will set it to something very low and tell you.
	std::cout << " StoppedMuonGun elow is too high " << elow_ << " resetting to 0 MeV" << std::endl;
      }
      
      spectrum_.initialize<MuonCaptureSpectrum>(elow_, ehi_,bin,kMaxUserSet,kMaxUser,kMaxMax,&randFlat_,&randUnitSphere_);

      double lowestEnergy = elow_;
      double upperEnergy  = ehi_;

      if (ehi_ > kMax) upperEnergy = kMax;
 
      double xLower = lowestEnergy/kMax;
      double xUpper = upperEnergy/kMax;
    
      //
      // integral of closure appoximation is 1/20 over [0,1].  this gives me the fraction of the spectrum we use, 
      // should weight overall rate by this

      fractionSpectrum_ = (20.) *  ( pow(xUpper,2)/2. - (4./3.)*pow(xUpper,3) + (7./4.)*pow(xUpper,4) - (6./5.)*pow(xUpper,5) 
			   + (1./3.)*pow(xUpper,6) )
	-  
	( pow(xLower,2)/2. - (4./3.)*pow(xLower,3) + (7./4.)*pow(xLower,4) - (6./5.)*pow(xLower,5) 
	  + (1./3.)*pow(xLower,6) );

      //
      // this is a DIFFERENT normalization.  Docdb 4378 and Armstrong et al tell us the rate about 57 MeV normalized to all
      // ordinary muon captures is 1.43 +-0.12 x 10^{-5}.  See the mathematica notebook in doc-db 16979.  Made configurable.

      omcNormalization_ = (rmcFrac)/(1/20. - 11432149083/pow(kMax,6) + (3610152342/5.)/pow(kMax,5) - (73892007/4.)/pow(kMax,4) + 246924/pow(kMax,3) - (3249/2.)/pow(kMax,2));
 
      if (physicsVerbosityLevel_ > 0){
	std::cout << "fraction of spectrum = " << fractionSpectrum_ << std::endl;
	std::cout << "rmc fraction         = " << rmcFrac           << std::endl;
	std::cout << "omc normalization    = " << omcNormalization_ << std::endl;
      }

    }
    else if (spectrumShape == "flat") {
      spectrum_.initialize<SimpleSpectrum>(elow_, ehi_, ehi_-elow_, SimpleSpectrum::Spectrum::Flat );
    }
    else {
      throw cet::exception("BADCONFIG")
        << "StoppedParticleMuonGun: unknown spectrum shape "<<spectrumShape<<"\n";
    }
  }

  //================================================================
  void StoppedMuonRMCGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    IO::StoppedParticleF stop = stops_.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);

    // next step is to get muon lifetime in the code together with the capture fraction 5/29/2018

    if (doHistograms_){
	_htZero->Fill(stop.t);
    }
    double energy = generateEnergy();

    double weight{0.};
    // two things can now happen with this photon.  It can proceed and possibly convert, or it can internally convert.
    // the probability of internal conversion is given by rho = 0.0069 (assume same as for RPC).  Throw a flat 

    if (randFlat_.fire() > rhoInternal_) {
      output->emplace_back( PDGCode::gamma, 
			    GenId::radiativeMuonCapture, 
			    pos,
			    CLHEP::HepLorentzVector(randUnitSphere_.fire(energy),energy), 
			    stop.t );

      event.put(std::move(output));

      // for future normalization

      double weightExternal = fractionSpectrum_ * omcNormalization_;
      std::unique_ptr<EventWeight> pw(new EventWeight(weightExternal));
      event.put(std::move(pw));
      weight = weightExternal;
      if ( doHistograms_ ) {
	_hmomentum->Fill(energy);
      }
    } 
    else {    // internal conversions

      CLHEP::HepLorentzVector mome, momp;

      muonCaptureSpectrum_.getElecPosiVectors(energy,mome,momp); 

      // CLHEP::HepLorentzVector fakeElectron( 105.*TMath::Sin(CLHEP::pi*60./180.),0.,105.*TMath::Cos(CLHEP::pi*60./180.),sqrt(105*105+me_*me_));
      // CLHEP::HepLorentzVector fakePositron(-105.*TMath::Sin(CLHEP::pi*60./180.),0.,105.*TMath::Cos(CLHEP::pi*60./180.),sqrt(105*105+me_*me_));

      output->emplace_back( PDGCode::e_minus, 
			    GenId::radiativeMuonCaptureInternal, 
			    pos,
			    mome, 
			    //fakeElectron, 
			    //			    800. );
      			    stop.t );
      output->emplace_back( PDGCode::e_plus, 
			    GenId::radiativeMuonCaptureInternal, 
			    pos,
			    momp,
			    //fakePositron, 
			    //			    800.);
      			    stop.t );

      event.put(std::move(output));

      // for future normalization
      double weightInternal = omcNormalization_*fractionSpectrum_*rhoInternal_;
      std::unique_ptr<EventWeight> pw(new EventWeight(weightInternal));
      weight = weightInternal;
      event.put(std::move(pw));
      if (verbosityLevel_ > 0) {
	std::cout << "original photon energy = " << energy << " and electron mass = " << me_ <<  std::endl;
	std::cout << "RMC electron/positron energies = " << mome.e() << " " << momp.e() << std::endl;
	std::cout << "and the full 4-vector: " << mome << " " << momp << std::endl;
	std::cout << "stop time = " << stop.t << std::endl;
	std::cout << " event weight = " << fractionSpectrum_ << " " << rhoInternal_ << " " << weightInternal << std::endl;
      }

      if ( doHistograms_ ) {
	_hWeight->Fill(weight);
	_hmomentum->Fill(energy);
	_hEnergyElectron->Fill(mome.e());
	_hEnergyPositron->Fill(momp.e());

	double mee = (mome+momp).m();
	_hMee->Fill(mee);
	_hMeeVsE->Fill(energy,mee);

	CLHEP::Hep3Vector p = mome.vect()+momp.vect();
	double y = (mome.e()-momp.e())/p.mag();

	_hy->Fill(y);

	// if ((mee > 85) && (energy > 85)) {
	//   printf(" >>> StoppedMuonRMCGun::produce EVENT: %10i energy: %12.5f mee: %12.5f\n",event.event(),energy,mee);
	// }

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
