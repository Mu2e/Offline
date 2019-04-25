//
// read ascii file for momentum and position of particle.  Set pdgid in constructor.
//
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

// random number includes
#include "SeedService/inc/SeedService.hh"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "CLHEP/Random/RandGaussQ.h"

// Framework includes
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"


// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "art_root_io/TFileService.h"


#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"

// C++ includes
#include <iostream>
#include <fstream>

using namespace std;

namespace mu2e {

  class FromAsciiMomentumAndPosition : public art::EDProducer {
  public:
    struct Config {
      fhicl::Atom<int> particlePdgId{ fhicl::Name("particlePdgId") };
      fhicl::Atom<std::string> vertexFileName{ fhicl::Name("vertexFileName") };
      fhicl::Atom<int> verbosityLevel{ fhicl::Name("verbosityLevel") };
      fhicl::Atom<bool> doHistograms{fhicl::Name("doHistograms"), false};
    };

    using Parameters = art::EDProducer::Table<Config>;

  private:
    int particlePdgId_;
    std::string vertexFileName_;
    int verbosityLevel_;
    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat    randFlat_;
    bool doHistograms_;

    ifstream vertexFile_;
    int lengthOfVertexFile_;
    int ncalls;



    double mass_;
    double protonMass_;

    TH1D* hIndexOfParticle;
    TH1D* hTimeOfParticle;
    TH1D* hRadiusOfParticle;
    TH1D* hMomentumOfParticle;

    std::vector<CLHEP::Hep3Vector> startPos_;
    std::vector<CLHEP::Hep3Vector> startMom_;
    std::vector<double> startTime_;
    std::vector<CLHEP::Hep3Vector> startInitialProtonMomentum_;
    std::vector<CLHEP::Hep3Vector> startInitialAntiProtonMomentum_;

    void produce(art::Event& event) override;

  public:
    explicit FromAsciiMomentumAndPosition(Parameters const& config);
  };

  FromAsciiMomentumAndPosition::FromAsciiMomentumAndPosition(Parameters const& config):
    EDProducer{config},
    //
    // Information from config file.
    particlePdgId_{config().particlePdgId()},
    vertexFileName_{config().vertexFileName()},
    verbosityLevel_{config().verbosityLevel()},
    eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())},
    randFlat_{eng_},
    doHistograms_{config().doHistograms()}
  {

    //pick up particle mass and other constants
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& p_data = pdt->particle(particlePdgId_).ref();
    mass_ = p_data.mass().value();

    const HepPDT::ParticleData& proton_data = pdt->particle(PDGCode::proton).ref();
    protonMass_ = proton_data.mass().value();

    ConfigFileLookupPolicy findConfig;
    std::string vertexFileString_ = findConfig(vertexFileName_);

   vertexFile_.open( vertexFileString_.c_str() );
   //
   // read them all in
   float xStart_,yStart_,zStart_,pxStart_,pyStart_,pzStart_, time_;
   float initialProtonMomentumX_,initialProtonMomentumY_,initialProtonMomentumZ_;
   float initialAntiProtonMomentumX_,initialAntiProtonMomentumY_,initialAntiProtonMomentumZ_;
   std::string line;

      while (std::getline(vertexFile_,line)){
       std::istringstream is(line);
       is >> xStart_ >> yStart_ >> zStart_ >> pxStart_ >> pyStart_ >> pzStart_ >> time_ >> initialAntiProtonMomentumX_ >>  initialAntiProtonMomentumY_ >> initialAntiProtonMomentumZ_
                   >> initialProtonMomentumX_ >>  initialProtonMomentumY_ >> initialProtonMomentumZ_;

       std::cout  << xStart_ << " " << yStart_ << " " << zStart_ << " " << pxStart_ << " " << pyStart_ << " " << pzStart_ << " " << time_
                   << " " << initialAntiProtonMomentumX_ << " " <<  initialAntiProtonMomentumY_ << " " << initialAntiProtonMomentumZ_
                  << " " << initialProtonMomentumX_ << " " <<  initialProtonMomentumY_ << " " << initialProtonMomentumZ_ << std::endl;
       startPos_.push_back(CLHEP::Hep3Vector(xStart_,yStart_,zStart_));
       startMom_.push_back(CLHEP::Hep3Vector(pxStart_,pyStart_,pzStart_));
       startTime_.push_back(time_);
       startInitialProtonMomentum_.push_back(
                                             CLHEP::Hep3Vector(initialProtonMomentumX_,initialProtonMomentumY_,initialProtonMomentumZ_));
       startInitialAntiProtonMomentum_.push_back(
                                                 CLHEP::Hep3Vector(initialAntiProtonMomentumX_,initialAntiProtonMomentumY_,initialAntiProtonMomentumZ_));
   }




   lengthOfVertexFile_ = startPos_.size();

      std::cout << "length of vertex file =  " << startPos_.size() << std::endl;

  ncalls = 0;


   produces <mu2e::GenParticleCollection>();

   if (verbosityLevel_ > 0){
     std::cout<<"FromAsciiMomentumAndPosition: Got particle pdgId and mass = " << particlePdgId_ << " "  << mass_<< ", input filename "<<config().vertexFileName()<<std::endl;
   }
    // set up histos
    if (doHistograms_)
      {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tfdir = tfs->mkdir( "FromAsciiMomentumAndPosition" );
        hIndexOfParticle = tfdir.make<TH1D>("hIndexOfParticle","index of chosen particle",2000,0.,2000.);
        hTimeOfParticle = tfdir.make<TH1D> ("hTimeOfParticle","particle time folded to accelerator cycle",200,0.,2000.);
        hRadiusOfParticle = tfdir.make<TH1D> ("hRadiusOfParticle","particle radius at starting z",100,0.,250.);
        hMomentumOfParticle = tfdir.make<TH1D> ("hMomentumOfParticle","particle momentum (mag) at starting z",2000,0.,2000.);
      }

  }

  void FromAsciiMomentumAndPosition::produce(art::Event& event )
  {

    auto output = std::make_unique<GenParticleCollection>();

     ++ncalls;

     //
     // pull a random entry from vertex file
     int ithPart = randFlat_.fire()*lengthOfVertexFile_;

     CLHEP::Hep3Vector momInitialParticle(startMom_.at(ithPart));
     CLHEP::Hep3Vector posInitialParticle(startPos_.at(ithPart));
     double eInitialParticle =
       sqrt( mass_*mass_ + momInitialParticle.mag()*momInitialParticle.mag() );
     CLHEP::HepLorentzVector fourMomInitialParticle
       (eInitialParticle,momInitialParticle);
     double timeTrial_ = startTime_.at(ithPart);
     CLHEP::Hep3Vector momInitialProton(startInitialProtonMomentum_.at(ithPart));
     double eInitialProton = sqrt(momInitialProton.mag()*momInitialProton.mag() + protonMass_*protonMass_);
     CLHEP::HepLorentzVector fourMomInitialProton(eInitialProton,momInitialProton);

     CLHEP::Hep3Vector momInitialAntiProton(startInitialAntiProtonMomentum_.at(ithPart));
     double eInitialAntiProton = sqrt(momInitialAntiProton.mag()*momInitialAntiProton.mag() + protonMass_*protonMass_);
     CLHEP::HepLorentzVector fourMomInitialAntiProton(eInitialAntiProton,momInitialAntiProton);

     if (timeTrial_ > 1695.){
       timeTrial_ = timeTrial_ - static_cast<int>(timeTrial_/1695.)*1695.;
     }
     if (doHistograms_){
       hIndexOfParticle->Fill(ithPart);
       hTimeOfParticle->Fill(timeTrial_);
       hRadiusOfParticle->Fill(sqrt((posInitialParticle.x()+3904.)*(posInitialParticle.x()+3904.) + posInitialParticle.y()*posInitialParticle.y()));
       hMomentumOfParticle->Fill(momInitialParticle.mag());
     }
    //
    // and now the particle itself
     output->emplace_back(PDGCode::type(particlePdgId_), GenId::fromAscii,
                          posInitialParticle,
                          fourMomInitialParticle,
                          timeTrial_);
     //
     // initial proton; doesn't matter where I put it, just want momentum vector
     output->emplace_back(PDGCode::type(PDGCode::proton), GenId::fromAscii,
                          CLHEP::Hep3Vector(0.,0.,0.),
                          fourMomInitialProton,
                          0. );
     //
     // initial antiproton; doesn't matter where I put it, just want momentum vector
     output->emplace_back(PDGCode::type(PDGCode::anti_proton), GenId::fromAscii,
                          CLHEP::Hep3Vector(0.,0.,0.),
                          fourMomInitialAntiProton,
                          0. );

    event.put(std::move(output));
    if (verbosityLevel_ > 0){
      std::cout << "looking at vertex file entry " << ithPart << std::endl;
      std::cout<< "FromAsciiMomentumAndPosition:  pdgid, position, momentum, time:" << "\n"
               << particlePdgId_ << "\n"
               << posInitialParticle << "\n"
               << fourMomInitialParticle << "\n"
               << timeTrial_ << std::endl;
      std::cout<< "FromAsciiMomentumAndPosition:  initialProton:" << "\n"
               << fourMomInitialProton << std::endl;
      std::cout<< "FromAsciiMomentumAndPosition:  initialAntiProton:" << "\n"
               << fourMomInitialAntiProton << std::endl;
    }
  }
  // FromAsciiMomentumAndPosition::generate

} //end namespace mu2e

DEFINE_ART_MODULE(mu2e::FromAsciiMomentumAndPosition);
