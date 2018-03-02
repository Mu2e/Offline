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
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
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
#include "art/Framework/Services/Optional/TFileService.h"


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

    TH1D* hIndexOfParticle;
    TH1D* hTimeOfParticle;
    TH1D* hRadiusOfParticle;
    TH1D* hMomentumOfParticle;
 
    std::vector<CLHEP::Hep3Vector> startPos_;
    std::vector<CLHEP::Hep3Vector> startMom_;
    std::vector<double> startTime_;
 
  public: 
    explicit FromAsciiMomentumAndPosition(Parameters const& config );
    virtual void produce(art::Event& event) override;
    virtual void beginRun(art::Run& ) override;
  };

  void FromAsciiMomentumAndPosition::beginRun(art::Run& ) {}

  FromAsciiMomentumAndPosition::FromAsciiMomentumAndPosition(Parameters const& config):
    //
    // Information from config file.
    particlePdgId_(config().particlePdgId()),
    vertexFileName_(config().vertexFileName()),
    verbosityLevel_(config().verbosityLevel()),
    eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    randFlat_(eng_),
    doHistograms_(config().doHistograms())
  {

    //pick up particle mass and other constants
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& p_data = pdt->particle(particlePdgId_).ref();
    mass_ = p_data.mass().value();


    ConfigFileLookupPolicy findConfig;
    std::string vertexFileString_ = findConfig(vertexFileName_);

   vertexFile_.open( vertexFileString_.c_str() );
   //
   // read them all in
   double xStart_,yStart_,zStart_,pxStart_,pyStart_,pzStart_, time_;

   do
     {
       vertexFile_ >> xStart_ >> yStart_ >> zStart_ >> pxStart_ >> pyStart_ >> pzStart_ >> time_;
       startPos_.push_back(CLHEP::Hep3Vector(xStart_,yStart_,zStart_));
       startMom_.push_back(CLHEP::Hep3Vector(pxStart_,pyStart_,pzStart_));
       startTime_.push_back(time_);
     } while (!vertexFile_.eof());

   lengthOfVertexFile_ = startPos_.size();
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

  //  FromAsciiMomentumAndPosition::~FromAsciiMomentumAndPosition(){}
  void FromAsciiMomentumAndPosition::produce(art::Event& event )
  {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

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
     output->push_back( GenParticle( PDGCode::type(particlePdgId_), GenId::fromAscii,
				    posInitialParticle,
				    fourMomInitialParticle,
				    timeTrial_ ));
   
    event.put(std::move(output));
    if (verbosityLevel_ > 0){
      std::cout<< "FromAsciiMomentumAndPosition:  pdgid, position, momentum, time:" << "\n"
	       << particlePdgId_ << "\n" 
	       << posInitialParticle << "\n"
	       << fourMomInitialParticle << "\n" 
	       << timeTrial_ << std::endl;
    }
  }
  // FromAsciiMomentumAndPosition::generate

} //end namespace mu2e

DEFINE_ART_MODULE(mu2e::FromAsciiMomentumAndPosition);
