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

#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

// ROOT includes
#include "TH1D.h"

// C++ includes
#include <iostream>

using namespace std;

namespace mu2e {

  class PrimaryAntiProtonGun : public art::EDProducer {

    fhicl::ParameterSet _pbarPhys;
    //
    // scheme for generation of pbars
    std::string _generatorType;

    //
    // which model we use
    std::string _generatorInvariantCrossSection;



    //
    // Offset of production point relative to the origin described above; in mm.
    CLHEP::Hep3Vector _beamDisplacementOnTarget;

    //
    // Rotation of beam direction wrt to target angle; in deg.
    double _beamRotationTheta;
    double _beamRotationPhi;
    double _beamRotationPsi;

    //
    // Beamspot is a 2D gaussian with this sigma in both x and y.
    double _beamSpotSigma;
  
    //
    // beam shape, flat or gaussian
    std::string _shape;
    double _rmax;


    CLHEP::HepRotation _gunRotation; 
    CLHEP::Hep3Vector _gunOrigin;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;
    double _elowFlat;
    double _ehiFlat;
 

    //
    // this is the most important one, probably the output of a separate G4 job that has protons interacting in the 
    // target
    std::string _inputProtonFile;

    //
    // where the proton hit points are if you want to fire from just one spot
    CLHEP::Hep3Vector _inputProtonPositionRelToTarget;



  
 

    bool _timeFolding;


    //
    // basic numbers for kinematics
    double _beamEnergy;
    double _beamMomentum;
    double _beamKineticEnergy;
 

    double _spectrumResolution;


    // Histogram control.
    bool _doHistograms;
    // Diagnostic histograms.
    TH1D* _hMultiplicity;
    TH1D* _hEPbar; 
    TH1D* _hCosThetaPbar;
    TH1D* _hradius;
    TH1D* _hzPosition;

    //
    // random numbers
    art::RandomNumberGenerator::base_engine_t& _eng;
    CLHEP::RandFlat    _randFlat;
    RandomUnitSphere   _randomUnitSphere;
    CLHEP::RandGaussQ  _randGaussQ;
 

    // Diagnostics
    int _diagLevel;

    //
    // end of configurable parameters

    BinnedSpectrum _pbarSpectrum;

    //
    // rotates target relative to Mu2e frame and the offset
    GenId::enum_type _pbarGenId;
    double _mass;

  public: 
    explicit PrimaryAntiProtonGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  PrimaryAntiProtonGun::PrimaryAntiProtonGun(const fhicl::ParameterSet& pset):
    //
    // Information from config file.
    _pbarPhys(pset.get<fhicl::ParameterSet>("physics")),
    _generatorType(_pbarPhys.get<string>("primaryantiprotonGun.type","flat")),
    _generatorInvariantCrossSection(_pbarPhys.get<string>("primaryantiprotonGun.invariantCrossSection","straightStriganov")),
    _beamDisplacementOnTarget(_pbarPhys.get<CLHEP::Hep3Vector>("beamDisplacementOnTarget")),
    _beamRotationTheta(_pbarPhys.get<double>("beamRotationTheta", 0)),
    _beamRotationPhi(_pbarPhys.get<double>("beamRotationPhi", 0)),
    _beamRotationPsi(_pbarPhys.get<double>("beamRotationPsi", 0)),
    _beamSpotSigma(_pbarPhys.get<double>("primaryProtonGun.beamSpotSigma")),    
    _shape(_pbarPhys.get<string>("primaryProtonGun.shape", "gaus")),
    _rmax(_pbarPhys.get<double>("primaryProtonGun.rmax", 100.)),

    _gunRotation(GeomHandle<ProductionTarget>()->protonBeamRotation()),
    _gunOrigin(GeomHandle<ProductionTarget>()->position()
               + _gunRotation*CLHEP::Hep3Vector(0., 0., GeomHandle<ProductionTarget>()->halfLength())),


    //
    // these configuration parameters go with a flat spectrum and assigns a weight
    _czmin(_pbarPhys.get<double>("primaryantiprotonGunFlat.cosmin",   -1.0)),
    _czmax(_pbarPhys.get<double>("primaryantiprotonGunFlat.cosmax",   1.0)),
    _phimin(_pbarPhys.get<double>("primaryantiprotonGunFlat.phimin",  -M_PI)),
    _phimax(_pbarPhys.get<double>("primaryantiprotonGunFlat.phimax",  +M_PI)),
    _elowFlat(_pbarPhys.get<double>("primaryantiprotonGunFlat.elow",      0.)),
    _ehiFlat(_pbarPhys.get<double>("primaryantiprotonGunFlat.ehi",        0.)),
 
    //
    // these configuration parameters go with a scheme where one throws flat and assigns a weight
    // or one chooses from the (momentum, angle) histograms 
    _inputProtonFile(_pbarPhys.get<string>("AntiProtonStartGenerator.Pbarfilename",
				      "ExampleDataFiles/ProtonInteractionsInTarget/protonInteractionsInTargetFromG4_00.txt")),

    _inputProtonPositionRelToTarget(_pbarPhys.get<CLHEP::Hep3Vector>("inputProtonPositionRelToTarget",CLHEP::Hep3Vector(0.,0.,0.))),
  
    _timeFolding(_pbarPhys.get<bool>("AntiProtonStartGenerator.foldingTimeOption", true)),
    _spectrumResolution(_pbarPhys.get<double>("spectrumResolution",0.)),
    //
    // everyone
    _doHistograms(_pbarPhys.get<bool>("primaryantiprotonGun.doHistograms", true)),
    // Histograms.
    _hEPbar(0),
    _hCosThetaPbar(0),
    _hradius(0),
    _hzPosition(0),

    //
    // random numbers
    _eng(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randFlat(_eng),    
    _randomUnitSphere(_eng, _czmin, _czmax, _phimin, _phimax ),
    _randGaussQ(_eng),
  

    // diagnostics
    _diagLevel (_pbarPhys.get<int>("primaryantiprotonGun.diagLevel",0)) {

    //pick up particle mass and other constants
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& p_data = pdt->particle(PDGCode::anti_proton).ref();
    _mass = p_data.mass().value();
    _beamMomentum = GlobalConstantsHandle<PhysicsParams>()->getProtonMomentum();
    _beamEnergy = GlobalConstantsHandle<PhysicsParams>()->getProtonEnergy();
    _beamKineticEnergy = GlobalConstantsHandle<PhysicsParams>()->getProtonKE();

    produces <mu2e::GenParticleCollection>();
 
    // Make ROOT subdirectory to hold diagnostic histograms; book those histograms.

    if ( _doHistograms ){
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "PrimaryAntiProton" );
      _hEPbar                  = tfdir.make<TH1D>( "hEElec",                  "PBAR Energy",                     100,    0.,   10.  );
      _hCosThetaPbar           = tfdir.make<TH1D>( "hCosThetaPbar",           "PBAR cos(theta)",                 100,    -1.,     1.   );
      _hradius       = tfdir.make<TH1D>( "hradius",       "PBAR radius from production target axis",         150,    0.,   150. );
      _hzPosition    = tfdir.make<TH1D>( "hzPosition",    "PBAR z Position from  production target center",  100, -100.,   100. );
    }

 
    // Pick correct model
    if ( _generatorInvariantCrossSection == "flat")
      { 
      _pbarGenId = GenId::pbarFlat;          
      _pbarSpectrum.initialize< SimpleSpectrum >( _elowFlat, _ehiFlat, _spectrumResolution, SimpleSpectrum::Spectrum::FlatTrunc );
    }
    else {
      throw cet::exception("MODEL")
        << "Wrong or not allowed PBAR energy spectrum";
    }
    
  }

  //  PrimaryAntiProtonGun::~PrimaryAntiProtonGun(){}

  void PrimaryAntiProtonGun::produce(art::Event& event ){

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    //
    //only allow one pbar; at 8 GeV primary, can't make two 
    //Loop over particles to generate

    //
    // these are needed for all schemes
    double time(0.);
    CLHEP::HepLorentzVector mom;
    CLHEP::Hep3Vector pos;
    CLHEP::HepRotation rot;
    double e(0.);
    double p(0.);
     if (_generatorType == "flat")
      {
 
	//
	// make these user-set so can tune in if needed.
	e = _elowFlat + _randFlat.fire();
	p = safeSqrt(e*e - _mass*_mass);


	// Simulate the size of the beam spot.
	double dx = 0;
	double dy = 0;
	double dr = 0;
	double phi = 0;

	if (_shape == std::string("gaus")) {
	  dx = _randGaussQ.fire();
	  dy = _randGaussQ.fire();
	}
	else if (_shape == std::string("flat")) {
	  // even in the circle
	  dr = _rmax*sqrt( _randFlat.fire() ); 
	  phi = 2.*M_PI*_randFlat.fire(); 
	  dx = dr*cos(phi);
	  dy = dr*sin(phi);
	}

	// Generated position.
	pos.set( _beamDisplacementOnTarget.x() + dx,
	         _beamDisplacementOnTarget.y() + dy,
	         _beamDisplacementOnTarget.z() );
	//
	// incoming angle
	rot.setTheta( _beamRotationTheta * CLHEP::degree );
	rot.setPhi( _beamRotationPhi * CLHEP::degree );
	rot.setPsi( _beamRotationPsi * CLHEP::degree );
	mom = rot * mom;
	//
	// Generated 4 momentum.
	CLHEP::HepLorentzVector mom(  _randomUnitSphere.fire(p), e );


	//
	// now we can put this antiproton in Mu2e coordinates
 
 
	if (_diagLevel > 0) 
	  {
	    std::cout << "from PrimaryAntiProtonGun, flat option, e and p = " << e << " " << p << std::endl;
	  }
      }

    else {
      throw cet::exception("MODEL")
        << "Wrong or unimplemented antiproton differential cross-section";
    }
     // Add the antiproton to the list of generated particles.
     // Convert position and momentum to Mu2e coordinates
     output->push_back( GenParticle( PDGCode::anti_proton, _pbarGenId,
				     _gunRotation*pos + _gunOrigin,
				     _gunRotation*mom,
				     time ));

     event.put(std::move(output));

    if( _doHistograms ){
      _hEPbar     ->Fill( e );
      _hCosThetaPbar        ->Fill( mom.cosTheta() );
      _hradius             ->Fill(sqrt( pos.x()*pos.x() + pos.y()*pos.y() ) );
      _hzPosition          ->Fill( pos.z() );
      }

  } // PrimaryAntiProtonGun::generate

} //end namespace mu2e

DEFINE_ART_MODULE(mu2e::PrimaryAntiProtonGun);
