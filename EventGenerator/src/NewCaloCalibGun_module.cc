// Sophie Middleton, 2022
// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <memory>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"

//For primary:
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//ROOT Includes
#include <TMath.h>
#include <TTree.h>
const double piconst = Pi();
using namespace std;
using namespace TMath;

namespace mu2e {

  //================================================================
  class NewCaloCalibGun : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> mean{Name("Mean"),1.},
      fhicl::Atom<double> energy{Name("PhotonEnergy"),6.13};//MeV
      fhicl::Atom<double> cosmin{Name("cosmin"),-1.},
      fhicl::Atom<double> cosmin{Name("cosmin"),1.},
      fhicl::Atom<double> phimin{Name("phimin"),0.},
      fhicl::Atom<double> phimin{Name("phimax"), CLHEP::twopi },
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), true },
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit NewCaloCalibGun(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:

    
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    GlobalConstantsHandle<ParticleDataTable> pdt;

    // Default values for the start and end of the live window.
    double _tmin = 0.;
    double _tmax = accPar->deBuncherPeriod;

    double _mean;
    double _energy;
    double _cosmin;
    double _cosmax;
    double _phimin;
    double _phimax;
    
    art::RandomNumberGenerator::base_engine_t& engine;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;

    const DetectorSystem  *_detSys;
    const DiskCalorimeter *_cal;
    int                    _nPipes;
    double                 _pipeRadius;
    std::vector<double>    _pipeTorRadius;
    std::vector<double>    _randomRad;
    CLHEP::Hep3Vector      _zPipeCenter;


    bool _doHistograms;
    TTree* _Ntupe;
    float _genErg;
    float _genTime;
    float _genCos;
    float _genPosX;
    float _genPosY;
    float _genPosZ;
 
  };

  //================================================================
  NewCaloCalibGun::NewCaloCalibGun(const Parameters& conf)
    : EDProducer{conf}
    , _mean(config.getDouble("caloCalibGun.mean",1.)),
    , _energy{conf().energy()}
    , _cosmin{conf().cosmin()}
    , _cosmax{conf().cosmax()}
    , _phimin{conf().phimin()}
    , _phimax{conf().phimax()}
    , _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , _randFlat{_engine},
    , _randPoissonQ{_engine, std::abs(_mean)},
    , _randomUnitSphere{_engine, _cosmin, _cosmax, 0, CLHEP::twopi},

  {
    produces<mu2e::StageParticleCollection>();
    produces<mu2e::PrimaryParticle>();
    
    _detSys          = &*GeomHandle<DetectorSystem>();
    _cal             = &*GeomHandle<DiskCalorimeter>();
    _nPipes          = _cal->caloInfo().getInt("nPipes");
    _pipeRadius      = _cal->caloInfo().getDouble("pipeRadius");
    _pipeTorRadius   = _cal->caloInfo().getVDouble("pipeTorRadius");
    _randomRad       = _cal->caloInfo().getVDouble("pipeTorRadius");
    //_zPipeCenter     = _cal->disk(0).geomInfo().origin()-CLHEP::Hep3Vector(0,0,_cal->disk(0).geomInfo().size().z()/2.0-_pipeRadius);
    // we normalize to the volume of the pipe (proportional to 2*pi*R if they have all the same radius) to draw a
    // random number from which to generate the photons
    double sumR(0);
    std::for_each(_randomRad.begin(), _randomRad.end(), [&](double& d) {sumR+=d; d = sumR; });
    std::for_each(_randomRad.begin(), _randomRad.end(), [&](double& d) {d /= sumR;});
    // Book histograms.
    if ( _doHistograms )
    {
      art::ServiceHandle<art::TFileService> tfs;
      _Ntupe  = tfs->make<TTree>("calibGun", "calibGun");
      _Ntupe -> Branch("genErg", &_genErg, "genErg/F");
      _Ntupe -> Branch("genTime", &_genTime, "genTime/F");
      _Ntupe -> Branch("genCos", &_genCos, "genCos/F");
      _Ntupe -> Branch("genPosX", &_genPosX, "genPosX/F");
      _Ntupe -> Branch("genPosY", &_genPosY, "genPosY/F");
      _Ntupe -> Branch("genPosZ", &_genPosZ, "genPosZ/F");
    }
    if (_mean < 0) throw cet::exception("RANGE") << "CaloCalibGun.mean must be non-negative "<< '\n';
  }

  //================================================================
  void NewCaloCalibGun::produce(art::Event& event) {
    auto output{std::make_unique<StageParticleCollection>()};
    //Call the gun here
    output->emplace_back(mustop,
                         ProcessCode::CaloCalib,
                         PDGCode::gamma,
                         mustop->endPosition(),
                         CLHEP::HepLorentzVector{randomUnitSphere_.fire(randomMom), randomE},
                         time
                         );
    
    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::NewCaloCalibGun);
