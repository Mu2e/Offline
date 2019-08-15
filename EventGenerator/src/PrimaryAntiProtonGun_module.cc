//
// Generate some number of antiprotons from proton collisions with target
// vertex only; antiprotons from daughters not permitted.
// based on algorithm in Striganov doc-db 1776, redone by Bernstein
// in doc-db 8448, etc.
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
#include "MCDataProducts/inc/ProcessCode.hh"


#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/copySimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"


// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <cstring>
#include <iostream>


#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "cetlib/map_vector.h"



using namespace std;

namespace mu2e {

  class PrimaryAntiProtonGun : public art::EDProducer {
  public:
    struct Config {
      fhicl::Atom<int> verbosityLevel{ fhicl::Name("verbosityLevel"),2 };
      fhicl::Atom<double> czmin{ fhicl::Name("czmin") , -1.};
      fhicl::Atom<double> czmax{ fhicl::Name("czmax") , 1.};
      fhicl::Atom<double> phimin{ fhicl::Name("phimin"), -M_PI};
      fhicl::Atom<double> phimax{ fhicl::Name("phimax"), +M_PI};
      fhicl::Atom<double> pLowFlat{ fhicl::Name("pLowFlat"), 0. };
      fhicl::Atom<double> pHighFlat{ fhicl::Name("pHighFlat"), 9. }; 
      fhicl::Atom<bool> doHistograms{fhicl::Name("doHistograms"), false};    
      fhicl::Atom<string> inputModuleLabel{fhicl::Name("inputModuleLabel"), "g4run"};    
    };

    using Parameters = art::EDProducer::Table<Config>;



  private:
    //
    // tmi line
    int _verbosityLevel;


    // Limits on angles
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    typedef cet::map_vector_key key_type;
 
    key_type                       idPbar;
    const unsigned                 stageOffsetPbar{0};
    const PDGCode::type            pdgIdPbar{PDGCode::anti_proton};
    const CLHEP::Hep3Vector&       positionPbar{0.,0.,0.};
    const CLHEP::HepLorentzVector& momentumPbar{0.,0.,0.,0.};
    double                         startGlobalTimePbar{0.};
    double                         startProperTimePbar{0.};
    unsigned                       startVolumeIndexPbar{0};
    unsigned                       startG4StatusPbar{1};

    // Limits on the generated energy.
    
    double _pLowFlat;
    double _pHighFlat;
    

    // random numbers
    art::RandomNumberGenerator::base_engine_t& _eng;
    CLHEP::RandFlat    _randFlat;
    RandomUnitSphere   _randomUnitSphere;
  
    double _spectrumResolution;
    BinnedSpectrum _pbarSpectrum;

    double _mass;

    bool _doHistograms;
    string _inputModuleLabel;

    TH1F* _hPbarMomentum;
    TH1F* _hPbarCosTheta;
    TH1F* _hPbarPhi;
    TH1F* _hNumberOfProtonInelastics;
    TH1F* _hNonProtonInelastics;

    TH1F* _hxLocationOfInelastic;
    TH1F* _hyLocationOfInelastic;
    TH1F* _hzLocationOfInelastic;
    TH1F* _hFoundInteractingProton;

    int numberOfProtonInelastics{0};
    int nonProtonInelastics{0};
    float eInitialProton{0.};
    bool first_protonInelastic{true};



  public: 
    explicit PrimaryAntiProtonGun(Parameters const& config );
    virtual void produce(art::Event& event) override;
    virtual void beginRun(art::Run& ) override;
  };

  void PrimaryAntiProtonGun::beginRun(art::Run& ) {}

  PrimaryAntiProtonGun::PrimaryAntiProtonGun(Parameters const& config):
    //
    // Information from config file.
    art::EDProducer{config},
    _verbosityLevel(config().verbosityLevel()),
    _czmin(config().czmin()),
    _czmax(config().czmax()),
    _phimin(config().phimin()),
    _phimax(config().phimax()),
    _pLowFlat(config().pLowFlat()),
    _pHighFlat(config().pHighFlat()),
   
    // these configuration parameters go with a scheme where one throws flat and assigns a weight
    // or one chooses from the (momentum, angle) histograms 
    // random numbers
    _eng(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randFlat(_eng),
    //
    // random numbers
    _randomUnitSphere(_eng, config().czmin(), config().czmax(), config().phimin(), config().phimax() ),
    _doHistograms(config().doHistograms()),
    _inputModuleLabel(config().inputModuleLabel())
  {

    //pick up particle mass and other constants
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& p_data = pdt->particle(PDGCode::anti_proton).ref();
    _mass = p_data.mass().value();



    produces <mu2e::SimParticleCollection>(""); 
    produces <mu2e::StepPointMCCollection>("");
 
    // set up histos
    if (_doHistograms)
      {
	art::ServiceHandle<art::TFileService> tfs;
	art::TFileDirectory tfdir = tfs->mkdir( "PrimaryAntiProtonGun" );

	_hPbarMomentum = tfdir.make<TH1F>("_hpBarMomentum","generated pbar momentum",100,0.,10000.);
	_hPbarCosTheta = tfdir.make<TH1F>("_hPbarCosTheta","antiproton cos theta in Lab Frame",100,-1.,1.);
	_hPbarPhi = tfdir.make<TH1F>("_hPbarPhi","antiproton phi in Lab Frame",100,-M_PI,+M_PI);
	_hNumberOfProtonInelastics = tfdir.make<TH1F>("_hNumberOfProtonInelastics","number of proton inelastics",10,0.,10.);
	_hNonProtonInelastics = tfdir.make<TH1F>("_hNonProtonInelastics","number of non-proton-inelastics but not killer volume",10,0.,10.);
	_hxLocationOfInelastic = tfdir.make<TH1F>("_hxLocationOfInelastic","x location of inelastic", 100,-100.,100.);
	_hyLocationOfInelastic = tfdir.make<TH1F>("_hyLocationOfInelastic","y location of inelastic", 100,-100.,100.);
	_hzLocationOfInelastic = tfdir.make<TH1F>("_hzLocationOfInelastic","z location of inelastic", 400,-6300,-5900.);//nominal is -6164.5+-80
	_hFoundInteractingProton = tfdir.make<TH1F>("_hFoundInteractingProton"," zero if wrote collection, one if did",2,0.,2.); 
      }

  }

  //  PrimaryAntiProtonGun::~PrimaryAntiProtonGun(){}
  void PrimaryAntiProtonGun::produce(art::Event& event )
  {


    //
    //only allow one pbar; at 8 GeV primary or less, can't make two with a 4.3 GeV threshold 

    //
    // declarations 
    CLHEP::Hep3Vector momInitialProton(0.,0.,0.);
    CLHEP::Hep3Vector posInitialProton(0.,0.,0.);
    CLHEP::HepLorentzVector momPbar(0.,0.,0.,0.);
    CLHEP::HepRotation rot;

    //
    // this routine looks through the simparticles in the event.  It finds the first proton inelastic collision
    // and writes out an antiproton thrown flat in 4pi wrt the momentum vector of the proton just before the inelastic
    // collision.  The momentum is chosen flat between 0 and the kinematic limit for the momentum of that proton.  

    art::Handle<SimParticleCollection> simParticleHandle;
    event.getByLabel(_inputModuleLabel,simParticleHandle);
    SimParticleCollection const& simParticles(*simParticleHandle);

    //
    // upgrade later for generality
    std::ostream& os = std::cout;

    if (_verbosityLevel > 0){
      os << "\n" << "SimParticleCollection has size = " << simParticles.size() << " particles." << std::endl;
      os << "ind      key    parent  pdgId       Start  Position            P            End Position               P     vol   process \n" ;
    }



    int numberOfProtonInelastics{0};
    int nonProtonInelastics{0};
    float eInitialProton{0.};
    bool first_protonInelastic{true};

    uint iInteracting{0};
    uint biggestSimParticleId{0};

    double xInelastic{-10000.};
    double yInelastic{-10000.};
    double zInelastic{0.};
 
    for (const auto& simPartPair: simParticles) {
      auto const& simPart = simPartPair.second;
      key_type const& key = simPartPair.first;

      if (simPart.id().asInt() > biggestSimParticleId){biggestSimParticleId = simPart.id().asInt();}

      art::Ptr<SimParticle> const& parentPtr = simPart.parent();
      int parentKey = -1;
      if(parentPtr) parentKey = int(parentPtr.key());
      if(_verbosityLevel > 0) {
	os
	  << std::setiosflags(std::ios::fixed | std::ios::right) 
	  << " " << std::setw(7) << key
	  << " " << std::setw(7) << parentKey
	  << " " << std::setw(8) << simPart.pdgId()
	  << " " << std::setw(8)  << std::setprecision(1) << simPart.startPosition().x()
	  << " " << std::setw(8)  << std::setprecision(1) << simPart.startPosition().y()
	  << " " << std::setw(8)  << std::setprecision(1) << simPart.startPosition().z()
	  << " " << std::setw(9) << std::setprecision(1) << simPart.startMomentum().vect().mag()
	  << "   "
	  << " " << std::setw(8)  << std::setprecision(1) << simPart.endPosition().x()
	  << " " << std::setw(8)  << std::setprecision(1) << simPart.endPosition().y()
	  << " " << std::setw(8)  << std::setprecision(1) << simPart.endPosition().z()
	  << " " << std::setw(9) << std::setprecision(1) << simPart.endMomentum().vect().mag()
	  << " " << std::setw(6) << simPart.endVolumeIndex()
	  << "  "
	  << " " << std::setiosflags(std::ios::left) << simPart.stoppingCode().name() 
	  << std::endl;
      }
      //
      // I want the proton at its first protonInelastic.  Eventually I will weight these events according to the 
      // measured total cross-section for pbar production so I can just take this. It could be wrong if the proton
      // has an inelastic collision and then survives for another inelastic collision that produces a pbar.
      // 
      // let's see how often this happens.  Note for some reason ProtonInelastic and protonInelastic are both process codes...
      // I don't think the former should ever happen so:
      if (simPart.stoppingCode() == ProcessCode::ProtonInelastic){
	throw cet::exception("ILLEGAL PROCESS CODE ProtonInelastic") << " from " << __func__ << std::endl;
      }
      if (simPart.stoppingCode() == ProcessCode::protonInelastic){ ++numberOfProtonInelastics;}
      if (simPart.stoppingCode() == ProcessCode::protonInelastic && first_protonInelastic){
	//
	// is key.asInt() by def'n = ind? seems that way 
	iInteracting = key.asInt(); 
	first_protonInelastic = false;
	eInitialProton = sqrt( simPart.endMomentum().vect().mag()*simPart.endMomentum().vect().mag() + _mass*_mass);
	xInelastic = simPart.endPosition().x();
	yInelastic = simPart.endPosition().y();
	zInelastic = simPart.endPosition().z();
	//
	//  need to save variables for populating pbar simParticle I will make later
	startVolumeIndexPbar =  simPart.endVolumeIndex();
      }

      if (simPart.stoppingCode() != ProcessCode::protonInelastic && simPart.stoppingCode() != ProcessCode::mu2eKillerVolume) {++nonProtonInelastics;}
    }

    //
    // and now 1/2 the real work: for the proton with the protonInelastic stopping code, create a new SimParticle with this proton information

    // first, we need to know the energy of the proton so we can set a limit for the energy of the antiproton
  

    _hNumberOfProtonInelastics->Fill(numberOfProtonInelastics);
    _hNonProtonInelastics->Fill(nonProtonInelastics);

    //
    // make new output collection; also return before I start defining kinematics on non-existent particle.  
    std::unique_ptr<mu2e::SimParticleCollection> outSimPartPtr(new mu2e::SimParticleCollection());
    //
    // will have to dereference the pointer to the simparticle collection
    auto& outSimPart = *outSimPartPtr;


    //
    // need a new StepPointMCCollection with the outgoing proton to trigger G4 in next stage
    std::unique_ptr<mu2e::StepPointMCCollection> outStepPointMCPtr(new mu2e::StepPointMCCollection());
    auto& outStepPoint = *outStepPointMCPtr;

    if(iInteracting <= 0) { // failed to find interacting proton
      //write empty collections
      event.put(std::move(outSimPartPtr), "");
      event.put(std::move(outStepPointMCPtr),"");
      _hFoundInteractingProton->Fill(0.5);
      return;	
    } else {_hFoundInteractingProton->Fill(1.5);}

    if (_verbosityLevel > 1) {
      std::cout << "interacting > 0" << xInelastic << " " << yInelastic << " " << zInelastic << std::endl;
    }

    //
    // need kinematic max for pbar in the Lab frame.
    // all the parameterizations we use assume a free proton at rest, even though the target is something else. Assume that for consistency.
    // 
    // the way I calculate the max energy in the Lab frame is:
    // 1) assume two protons smashing together along the z axis; one is our beam proton, which may have lost some energy before interacting, one the target at rest
    // 2) calculate the maximum pbar energy for that beam proton
    // 3) boost back to lab.  I only care about the magnitude of the maximum momentum for this calculation, so none of the directions matter!  
  
    // first, what is the maximum energy possible in the cm frame: need some kinematic variables
    double _mandelS = 2.*_mass*eInitialProton + _mass*_mass + _mass*_mass;
    double _rootS = sqrt(_mandelS);
    double protonEnergyInCM_2 = (_rootS/2.)*(_rootS/2.);

    //
    //these quantities belong to the antiproton; you can look these up (one of many places: F. Taylor et al., PRD14,No. 5, 1 Sept 1976 p. 1217)
    double eStarMaxInCM_ = (_mandelS + _mass*_mass - (3.*_mass)*(3.*_mass))/(2.*sqrt(_mandelS));
    double pStarMaxInCM_ = sqrt(eStarMaxInCM_*eStarMaxInCM_ - _mass*_mass);
    CLHEP::HepLorentzVector maxAntiProtonInCM(CLHEP::Hep3Vector(0.,0.,pStarMaxInCM_),eStarMaxInCM_);
 
    //
    // here is the initial incoming proton in the CM frame
    CLHEP::HepLorentzVector initialProtonInCM(0.,0.,sqrt(protonEnergyInCM_2 - _mass*_mass),_rootS/2.);
    CLHEP::Hep3Vector betaCM = initialProtonInCM.boostVector();

    //
    // take the biggest possible antiproton energy in the CM and boost it back to the Lab. betaCM sign is defined by CLHEP, ultimately same as ROOT.
    double pHighThisProton_ = maxAntiProtonInCM.boost(betaCM).vect().mag();



 
    if (_verbosityLevel > 1)
      {
	std::cout << "mandelstam S = " << _mandelS << std::endl;
	std::cout << "max energy in cm, max momentum in cm = " << eStarMaxInCM_ << " " << pStarMaxInCM_ << std::endl;
	std::cout << " initial proton energy and maximum pbar momentum = " <<  eInitialProton << " " << pHighThisProton_ << std::endl;
      } 
    double pPbar = _randFlat.fire(_pLowFlat,pHighThisProton_);
    double ePbar = sqrt(_mass*_mass + pPbar*pPbar);

    if (_verbosityLevel > 0)
      {
	std::cout << " pLowFlat = " << _pLowFlat << " upperMomentum " << pHighThisProton_ << " " << ePbar << std::endl;
      }

    //
    // Generated  momentum.
    momPbar.setVect(_randomUnitSphere.fire(pPbar));
    momPbar.setE(ePbar);

    if (_verbosityLevel > 0) 
      {
	std::cout << "from PrimaryAntiProtonGun, flat option, e and p = " << ePbar << " " << pPbar << " " << momPbar << std::endl;
      }
  

    // Add the antiproton to the list of generated particles.
    // Convert position and momentum to Mu2e coordinates
    if (_doHistograms)
      {
	_hPbarMomentum->Fill(momPbar.vect().mag());
	_hPbarCosTheta->Fill(momPbar.cosTheta());
	_hPbarPhi->Fill(momPbar.phi());
	_hxLocationOfInelastic->Fill(xInelastic-3904.);
	_hyLocationOfInelastic->Fill(yInelastic);
	_hzLocationOfInelastic->Fill(zInelastic);
	if (_verbosityLevel > 0)
	  {
	    std::cout << " momentum, cos Theta, phi = " << momPbar.vect().mag() << " " << momPbar.cosTheta() << " " << momPbar.phi() << std::endl;
	  }

      }


    //
    // in order to make the pbar's art pointer to parent, need to 
    // collect info about the parent product
    auto SPpid = event.getProductID<SimParticleCollection>("");
    //      auto SPpid = simParticleHandle.id();
    auto SPpg  = event.productGetter(SPpid);
    copySimParticleCollection(simParticles,outSimPart,SPpid,SPpg);

    //
    // make the parent art ptr
    art::Ptr<SimParticle> pptr(SPpid,size_t(iInteracting),SPpg);
    //
    // make new SimParticle:  create pbar, copy starting point..


    if (_verbosityLevel > 0){
      os << "size of sim particle collection before pbar = " << outSimPart.size() << std::endl;
    }

    // get key for new simpart
    // https://internal.dunescience.org/doxygen/map__vector_8h_source.html for next line; 
    key_type newPbarKey = cet::map_vector_key(biggestSimParticleId + 1);


    if (_verbosityLevel > 1){
      os << "made newPbarKey, about to make newPbar" << std::endl;
    }
    auto const& oldParent = simParticles[key_type(iInteracting)]; 
    SimParticle newPbar(
			newPbarKey                            // id
			,0                                     // stageOffset
			,pptr                                  // parentSim
			,PDGCode::anti_proton                  // pdgId
			,art::Ptr<GenParticle>()               // since this comes from a SimParticle the ptr to a GenParticle should be null
			,oldParent.endPosition()                   // born where parent stopped
			,momPbar                              // momentum created here
			,oldParent.endGlobalTime()                 // pbar created at same time as parent interacted
			,0.                                    // proper time is zero for this new particle
			,oldParent.endVolumeIndex()                // starts where parent ends
			,1                                     // G4 status for new particle
			,ProcessCode::mu2eProtonInelastic      // defined for this sort of process
			);
    if (_verbosityLevel > 1){
      os << "made newPbar" << std::endl;
    }


    //
    // add end information to complete SimParticle
    newPbar.addEndInfo(
		       oldParent.endPosition()                // hasn't gone anywhere
		       ,momPbar                               // ditto, just fill in
		       ,oldParent.endGlobalTime()             // ditto
		       ,0.                                    // again no change in proper time
		       ,oldParent.endVolumeIndex()            // same
		       ,1                                     // G4 status for where pbar ends
		       ,ProcessCode::mu2eProtonInelastic      // same process code
		       ,momPbar.e()- _mass                    // preLastStepKE
		       ,0                                     // nSteps
		       ,-1.                                   // track length
		       );                 
    //
    // next stage will fail unless you tell it this new particle is the daughter of the proton!
    art::Ptr<SimParticle> pbarAsDaughter(SPpid,size_t(newPbarKey.asInt()),SPpg);
    auto & newCopyParent = outSimPart[key_type(iInteracting)];
    newCopyParent.addDaughter(pbarAsDaughter);
 
    if (_verbosityLevel > 1){
      os << "made newPbar and added endInfo" << std::endl;
    }
 
    outSimPart[newPbarKey] = newPbar;


    //
    // write out new collection

    if (_verbosityLevel > 1){
      os << " \n \n \n" << "output simpart: " << "\n" << std::endl;

      for (const auto& outPartPair: outSimPart){
	auto const& outPart = outPartPair.second;
	key_type const outKey = outPartPair.first;
	art::Ptr<SimParticle> const& outParentPtr = outPart.parent();
	int outParentKey = -1;
	if(outParentPtr) outParentKey = int(outParentPtr.key());
	os 
	  << std::setiosflags(std::ios::fixed | std::ios::right) 
	  << " " << std::setw(7) << outKey
	  << " " << std::setw(7) << outParentKey
	  << " " << std::setw(8) << outPart.pdgId()
	  << " " << std::setw(8)  << std::setprecision(1) << outPart.startPosition().x()
	  << " " << std::setw(8)  << std::setprecision(1) << outPart.startPosition().y()
	  << " " << std::setw(8)  << std::setprecision(1) << outPart.startPosition().z()
	  << " " << std::setw(9) << std::setprecision(1) << outPart.startMomentum().vect().mag()
	  << "   "
	  << " " << std::setw(8)  << std::setprecision(1) << outPart.endPosition().x()
	  << " " << std::setw(8)  << std::setprecision(1) << outPart.endPosition().y()
	  << " " << std::setw(8)  << std::setprecision(1) << outPart.endPosition().z()
	  << " " << std::setw(9) << std::setprecision(1) << outPart.endMomentum().vect().mag()
	  << " " << std::setw(6) << outPart.endVolumeIndex()
	  << "  "
	  << " " << std::setiosflags(std::ios::left) << outPart.stoppingCode().name() 
	  << std::endl;
     
      }
    }

    //
    // I need to create a StepPointMC
    bool foundThePbar{false};
    for (const auto& outPartPair: outSimPart){
      auto const& outPart = outPartPair.second;
      key_type const outKey = outPartPair.first;
      if (_verbosityLevel > 1){
	os << "copied collection PDGid = " << outPart.pdgId() << std::endl;
      }
      if (outPart.pdgId() == PDGCode::anti_proton){
	if (foundThePbar){
	  throw cet::exception("FOUND MULTIPLE PBARS IN COLLECTION") << " from " << __func__ << std::endl;
	}
	foundThePbar = true;
	//
	// make the parent art ptr
	if (_verbosityLevel > 1){
	  os << " about to make pptr for newStepPointMC" << std::endl;
	}

	art::Ptr<SimParticle> pptr(SPpid,size_t(outKey.asInt()),SPpg);

	if (_verbosityLevel > 1){
	  os << " made pptr for newStepPointMC" << std::endl;
	}

	StepPointMC newStepPointMC(
				   pptr                                  // SimPart constructor with an art::Ptr
				   ,oldParent.endVolumeIndex()           // starts where parent ends
				   ,0.                                   // total energy deposit
				   ,0.                                   // non-ionizing energy deposit
				   ,oldParent.endGlobalTime()            // as above
				   ,0.                                   // proper time is zero for this new particle
				   ,oldParent.endPosition()              // born where parent stopped
				   ,momPbar.vect()                       // momentum created here
				   ,0.                                   // step length
				   ,ProcessCode::Transportation          // this is like a virtual detector
				   );
	outStepPoint.push_back(newStepPointMC);  //not a map just a regular std::vector...yay...
      }
    }
    if (!foundThePbar) {
      throw cet::exception("DIDNT FIND PBAR IN COLLECTION") << " from " << __func__ << std::endl;
    }


    event.put(std::move(outSimPartPtr), "");
    event.put(std::move(outStepPointMCPtr),"");
  }


// PrimaryAntiProtonGun::generate

} //end namespace mu2e

DEFINE_ART_MODULE(mu2e::PrimaryAntiProtonGun);
