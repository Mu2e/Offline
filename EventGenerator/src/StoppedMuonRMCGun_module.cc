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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Sequence.h"

// Mu2e includes
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/Mu2eUtilities/inc/Table.hh"
#include "Offline/Mu2eUtilities/inc/RootTreeSampler.hh"
#include "Offline/GeneralUtilities/inc/RSNTIO.hh"

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
    const double czmax_;
    const double czmin_;
    std::unique_ptr<CLHEP::RandGeneral>  randSpectrum_;
    RandomUnitSphere     randUnitSphere_;
    RandomUnitSphere     randUnitSphereExt_; //For photons, to limit cosz
    CLHEP::RandFlat      randFlat_;

    MuonCaptureSpectrum  muonCaptureSpectrum_;

    RootTreeSampler<IO::StoppedParticleF> stops_;

    bool   kMaxUserSet_;
    double kMaxUser_;
    bool   doHistograms_;

    double me_;                        // electron mass
    double mmu_;

    double fractionSpectrum_;
    double omcNormalization_;
    double internalNormalization{0.};
    double externalNormalization{0.};

    TH1F* _hmomentum;
    TH1F* _hCosz;
    TH1F* _hEnergyElectron;
    TH1F* _hEnergyPositron;
    TH1F* _hWeight;
    TH1F* _htZero;
    TH1F* _hMee;
    TH2F* _hMeeVsE;
    TH1F* _hMeeOverE;
    TH1F* _hy;                                // splitting function
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
    double generateEnergy();
    double integrateClosure(const double xLow, const double xHigh);

    void   parseSpectrumShape(const fhicl::ParameterSet& psphys);

  public:
    explicit StoppedMuonRMCGun(const fhicl::ParameterSet& pset);
    ~StoppedMuonRMCGun();
    virtual void produce(art::Event& event);
  };

//================================================================
  StoppedMuonRMCGun::StoppedMuonRMCGun(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , psphys_             (pset.get<fhicl::ParameterSet>("physics"))
    , rhoInternal_        (psphys_.get<double>("rhoInternal"))
    , spectrum_           (BinnedSpectrum(psphys_))
    , verbosityLevel_     (pset.get<int>("verbosityLevel", 0))
    , eng_                (createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , czmax_              (pset.get<double>("czmax",  1.))
    , czmin_              (pset.get<double>("czmin", -1.))
    , randUnitSphere_     (eng_)
    , randUnitSphereExt_  (eng_, czmax_, czmin_)
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

    me_  = GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::e_minus ).mass();
    mmu_ = GlobalConstantsHandle<ParticleDataList>()->particle(PDGCode::mu_minus).mass();

    // initialize binned spectrum - this needs to be done right
    parseSpectrumShape(psphys_);

    randSpectrum_ = std::make_unique<CLHEP::RandGeneral>(eng_, spectrum_.getPDF(), spectrum_.getNbins());

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "StoppedMuonRMCGun" );

      _hmomentum       = tfdir.make<TH1F>("hmomentum", "Produced photon momentum, RMC", 70,  0.,  140.  );
      _hCosz           = tfdir.make<TH1F>("hCosz", "Produced external photon cosz, RMC", 200,  -1.,  1.  );
      _hEnergyElectron = tfdir.make<TH1F>("hEnergyElectron", "Produced electron energy, RMC Internal", 70,  0.,  140.  );
      _hEnergyPositron = tfdir.make<TH1F>("hEnergyPositron", "Produced electron energy, RMC Internal", 70,  0.,  140.  );
      _htZero          = tfdir.make<TH1F>("htZero"         , "Stopped Muon time", 100,0.,2000.);
      _hWeight         = tfdir.make<TH1F>("hWeight"        , "Event Weight ", 100,0.,2.e-5);
      _hMee            = tfdir.make<TH1F>("hMee"           , "M(e+e-) "     , 200,0.,200.);
      _hMeeVsE         = tfdir.make<TH2F>("hMeeVsE"        , "M(e+e-) "     , 200,0.,200.,200,0,200);
      _hMeeOverE       = tfdir.make<TH1F>("hMeeOverE"      , "M(e+e-)/E"          , 200, 0.,1);
      _hy              = tfdir.make<TH1F>("hy"             , "y = (ee-ep)/|pe+pp|", 200,-1.,1.);
    }
  }

  //================================================================
  StoppedMuonRMCGun::~StoppedMuonRMCGun() {
  }

  //================================================================
  void StoppedMuonRMCGun::parseSpectrumShape(const fhicl::ParameterSet& psphys) {

    const std::string spectrumShape(psphys.get<std::string>("spectrumShape"));
    const int physicsVerbosityLevel_(psphys.get<int>("physicsVerbosityLevel"));

    if (spectrumShape == "RMC") {

      bool   blind       = psphys.get<bool>  ("blind");
      bool   kMaxUserSet = psphys.get<bool>  ("kMaxUserSet");

      double kMaxUser(0);
      if (kMaxUserSet) kMaxUser = psphys.get<double>("kMaxUser");

      double rmcFrac     = psphys.get<double>("rmcFrac");

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

      if ( spectrum_.getXMin() > kMax ) {
        //
        // if I told you what kMax was you could unblind kMax.  Therefore I will set it to something very low and tell you.
        std::cout << " StoppedMuonGun elow is too high " << spectrum_.getXMin() << " resetting to 0 MeV" << std::endl;
      }

      double lowestEnergy = spectrum_.getXMin();
      double upperEnergy  = spectrum_.getXMax();

      if (spectrum_.getXMax() > kMax) upperEnergy = kMax;
      // papers measure R(photon>57) = 1.43e-05. Hardwire that.
      const double rGammaEnergy = 57.; // this is what was measured, won't change unless someone does it again. Measurements are e>57.

      if (spectrum_.getXMin() < rGammaEnergy){
        lowestEnergy = rGammaEnergy;
        std::cout << "inside " << __func__ << " resetting lower energy to physical limit from " << spectrum_.getXMin() << " to " << rGammaEnergy << std::endl;
      }
      if (spectrum_.getXMax() > kMax) {
        upperEnergy = kMax;
        std::cout << "inside " << __func__ << " resetting upper energy to physical limit from " << spectrum_.getXMax() << " to " << kMax << std::endl;
      }

      double xLower = lowestEnergy/kMax;
      double xUpper = upperEnergy/kMax;
      const double xGammaEnergy = rGammaEnergy/kMax;

      //
      // closure approximation is R(photon>57 MeV) = ( e^2/pi)*(kMax/muonMass)^2*(1 - (N-Z)/(N+Z))* integral from 57/kmax to 1 of (1 -2x + 2x^2)x(1-x)^2 dx
      // and the Bergsbusch et al paper docdb 1192 says for Al this is measured to be 1.43 times 10^-5. Now the integral above varies with kmax.  I am going to pin the
      // integral to the data. So normalization is (integral from elow to ehi / integral from 57/kmax to 1), or the fraction of the area we look at, times 1.43 x 10^-5

      double fractionOfSpectrum = integrateClosure(xLower,xUpper)/integrateClosure(xGammaEnergy,1.);

      externalNormalization = fractionOfSpectrum*rmcFrac*(czmax_-czmin_)/2.; //including fraction of cosz simulated in externals
      internalNormalization = rhoInternal_*fractionOfSpectrum*rmcFrac;

      if (physicsVerbosityLevel_ > 0){
        std::cout << "lowestEnergy, upperEnergy, xLower, xUpper, xGammaEnergy, kMax, rmcFrac, externalNormalization, internalNormalization" << "\n" <<
          lowestEnergy<< " " << upperEnergy<< " " << xLower<< " " << xUpper<< " " << xGammaEnergy<< " "
                  << kMax<< " " << rmcFrac<< " " << externalNormalization<< " " << internalNormalization  << std::endl;
        std::cout << "fraction of spectrum = " << fractionOfSpectrum << std::endl;
      }
    }
  }

  //================================================================
  double StoppedMuonRMCGun::generateEnergy() {
    return spectrum_.sample(randSpectrum_->fire());
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
      CLHEP::HepLorentzVector photon(randUnitSphereExt_.fire(energy),energy);
      output->emplace_back( PDGCode::gamma,
                            GenId::ExternalRMC,
                            pos,
                            photon,
                            stop.t );

      event.put(std::move(output));
      std::unique_ptr<EventWeight> pw(new EventWeight(externalNormalization));
      event.put(std::move(pw));

      weight = externalNormalization;
      if ( doHistograms_ ) {
        _hmomentum->Fill(energy);
        _hCosz->Fill(photon.cosTheta());
        _hWeight->Fill(weight);
      }
    }
    else {    // internal conversions

      CLHEP::HepLorentzVector mome, momp;

      muonCaptureSpectrum_.getElecPosiVectors(energy,mome,momp);

      // CLHEP::HepLorentzVector fakeElectron( 105.*TMath::Sin(CLHEP::pi*60./180.),0.,105.*TMath::Cos(CLHEP::pi*60./180.),sqrt(105*105+me_*me_));
      // CLHEP::HepLorentzVector fakePositron(-105.*TMath::Sin(CLHEP::pi*60./180.),0.,105.*TMath::Cos(CLHEP::pi*60./180.),sqrt(105*105+me_*me_));

      output->emplace_back( PDGCode::e_minus,
                            GenId::InternalRMC,
                            pos,
                            mome,
                            //fakeElectron,
                            //                            800. );
                            stop.t );
      output->emplace_back( PDGCode::e_plus,
                            GenId::InternalRMC,
                            pos,
                            momp,
                            //fakePositron,
                            //                            800.);
                            stop.t );

      event.put(std::move(output));

      // for future normalization
      std::unique_ptr<EventWeight> pw(new EventWeight(internalNormalization));
      event.put(std::move(pw));

      if (verbosityLevel_ > 0) {
        std::cout << "original photon energy = " << energy << " and electron mass = " << me_ <<  std::endl;
        std::cout << "RMC electron/positron energies = " << mome.e() << " " << momp.e() << std::endl;
        std::cout << "and the full 4-vector: " << mome << " " << momp << std::endl;
        std::cout << "stop time = " << stop.t << std::endl;
        std::cout << " event weight = " << fractionSpectrum_ << " " << rhoInternal_ << " " << internalNormalization << std::endl;
      }
      weight = internalNormalization;
      if ( doHistograms_ ) {
        _hWeight->Fill(weight);
        _hmomentum->Fill(energy);
        _hEnergyElectron->Fill(mome.e());
        _hEnergyPositron->Fill(momp.e());

        double mee = (mome+momp).m();
        _hMee->Fill(mee);
        _hMeeVsE->Fill(energy,mee);
        _hMeeOverE->Fill(mee/energy);

        CLHEP::Hep3Vector p = mome.vect()+momp.vect();
        double y = (mome.e()-momp.e())/p.mag();

        _hy->Fill(y);
      }
    }
  }

  //================================================================
  double StoppedMuonRMCGun::integrateClosure(const double xLow, const double xHigh){
    const double xHi2 = xHigh*xHigh;
    const double xLow2 = xLow*xLow;
    double result = (xHi2)/2. - (4./3.)*xHi2*xHigh + (7./4.)*(xHi2*xHi2) - (6./5.)*(xHi2)*(xHi2)*xHigh + (1./3.)*(xHi2*xHi2*xHi2)
      - ( (xLow*xLow)/2. - (4./3.)*xLow2*xLow + (7./4.)*(xLow2*xLow2) - (6./5.)*(xLow2)*(xLow2)*xLow + (1./3.)*(xLow2*xLow2*xLow2) );
    return result;
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedMuonRMCGun)
