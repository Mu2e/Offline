// Anthony Palladino, 2016
//
// Derived from an earlier file from rhbob, now inherits from
// a different base class (EDProducer).
//
// Generate muonic-Aluminum X-rays from the muon stop distribution.

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
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

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
#include "CLHEP/Random/RandFlat.h"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;

namespace mu2e {

  //================================================================
  class StoppedMuonXRayGammaRayGun : public art::EDProducer {
    fhicl::ParameterSet _psphys;

    // MuonicXRay momentum.
    double _p;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    art::RandomNumberGenerator::base_engine_t& _eng;
    RandomUnitSphere _randomUnitSphere;
    CLHEP::RandFlat _randFlat;
    CLHEP::RandExponential  _randExp;
    RootTreeSampler<IO::StoppedParticleF> _stops;

    // Control which photons we want to simulate
    bool _do66;
    bool _do347;
    bool _do844;
    bool _do1809;

    // Control histograms.
    bool _doHistograms;

    // Histograms.
    TH1F* _hMultiplicity;
    TH1F* _hcz;
    TH1F* _hphi;
    TH1F* _hmomentum;
    TH1F* _hradius;
    //TH1F* _hzPos;
    TH1F* _htime;
    TH2F* _hxyPos;
    //TH2F* _hrzPos;

    void bookHistograms();

  public:
    explicit StoppedMuonXRayGammaRayGun(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  StoppedMuonXRayGammaRayGun::StoppedMuonXRayGammaRayGun(const fhicl::ParameterSet& pset):
    EDProducer{pset},
    _psphys(pset.get<fhicl::ParameterSet>("physics")),
    _p(0.0),
    _czmin (_psphys.get<double>("czmin", -1.0)),
    _czmax (_psphys.get<double>("czmax",  1.0)),
    _phimin(_psphys.get<double>("phimin", 0. )),
    _phimax(_psphys.get<double>("phimax", CLHEP::twopi )),
    _eng(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randomUnitSphere(_eng, _czmin, _czmax, _phimin, _phimax ),
    _randFlat(_eng),
    _randExp(_eng),
    _stops(_eng, pset.get<fhicl::ParameterSet>("muonStops")),
    _do66(_psphys.get<bool>("do66", true )),
    _do347(_psphys.get<bool>("do347", true )),
    _do844(_psphys.get<bool>("do844", true )),
    _do1809(_psphys.get<bool>("do1809", true )),
    _doHistograms(_psphys.get<bool>("doHistograms", true )),
    _hMultiplicity(0),
    _hcz(0),
    _hphi(0),
    _hmomentum(0),
    _hradius(0),
    //_hzPos(0),
    _htime(0),
    _hxyPos(0){ //,_hrzPos(0)

    produces<mu2e::GenParticleCollection>();

    if ( _doHistograms ) bookHistograms();
  }


  //================================================================
  void StoppedMuonXRayGammaRayGun::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const auto& stop = _stops.fire();

    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    double time = stop.t;
    const double genRadius = sqrt((stop.x+3904.0)*(stop.x+3904.0)+stop.y*stop.y);

    int nphotons = 0;
    vector<double> photonEnergy;
    vector<double> photonTime;

    // create X Rays and Gamma Rays:
    double prob = _randFlat.fire();
    if (_do66 && prob < 0.625){  // 3d-2p line
      ++nphotons;
      photonEnergy.push_back(0.0661);
      photonTime.push_back(time);
    }
    prob = _randFlat.fire();
    if (_do347 && prob < 0.798){  // 2p-1s line
      ++nphotons;
      photonEnergy.push_back(0.3468);
      photonTime.push_back(time);
    }
    prob = _randFlat.fire();
    if (_do844 && prob < 0.040){  //
      ++nphotons;
      photonEnergy.push_back(0.844);
      //Note: This is a delayed gamma, need to add delay time
      double meanLifetime844 = 822.0*CLHEP::second; //822 second lifetime (same as 9.5min(570s) halflife)
      photonTime.push_back(time + _randExp.fire(meanLifetime844));

    }
    prob = _randFlat.fire();
    if (_do1809 && prob < 0.300){  //
      ++nphotons;
      photonEnergy.push_back(1.809);
      //Note: This is a semi-prompt gamma, need to add delay time
      double meanLifetime1809 = 864.0*CLHEP::ns; //864ns, same lifetime as muonic Aluminum
      photonTime.push_back(time + _randExp.fire(meanLifetime1809));
    }

    for (int ithphoton=0; ithphoton < nphotons; ++ithphoton){
      // Compute momentum 3-vector
      CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(photonEnergy[ithphoton]);

      // Compute energy
      _p = photonEnergy[ithphoton];
      double e = _p; // yes this is stupid, let the optimizer fix it.  keeps code parallel among guns

      // Set four-momentum
      CLHEP::HepLorentzVector mom(p3, e);

      //time for this photon
      double timephoton = photonTime[ithphoton];

      // Add the particle(s) to  the list.
      output->push_back( GenParticle( PDGCode::gamma,
                                     GenId::StoppedMuonXRayGammaRayGun, pos, mom, timephoton));

      // Fill histograms
      if (_doHistograms){
          _hMultiplicity->Fill(nphotons);
          _hcz->Fill(p3.cosTheta());
          _hphi->Fill(p3.phi());
          _hmomentum->Fill(_p);
          _hradius->Fill( genRadius );
          //_hzPos->Fill(stop.z);
          _htime->Fill(timephoton);
          _hxyPos->Fill( stop.x+3904.0, stop.y   );
          //_hrzPos->Fill( stop.z, genRadius );
      }
    }


    event.put(std::move(output));
  }


  void StoppedMuonXRayGammaRayGun::bookHistograms(){

    // Compute a binning that ensures that the stopping target foils are at bin centers.
    //GeomHandle<StoppingTarget> target;
    //Binning bins = zBinningForFoils(*target,7);
    //Binning bins2 = zBinningForFoils(*target,3);

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir( "StoppedMuonXRayGammaRayGun" );

    _hMultiplicity = tfdir.make<TH1F>( "hMultiplicity",
                                       "MuonicXRay Multiplicity",
                                       10,  0.,  10.  );
    _hcz           = tfdir.make<TH1F>( "hcz",
                                       "MuonicXRay Photon cos(theta) at Production;(MeV)",
                                       100,  -2.,  2.  );
    _hphi          = tfdir.make<TH1F>( "hphi",
                                       "MuonicXRay Photon phi at Production;(MeV)",
                                       100,  -M_PI,  M_PI  );
    _hmomentum     = tfdir.make<TH1F>( "hmomentum",
                                       "MuonicXRay Photon Momentum at Production;(MeV)",
                                       1000,  0.,  2.0  );
    _hradius       = tfdir.make<TH1F>( "hradius",
                                       "MuonicXRay Photon Radius at Production;(mm)",
                                       60,  0., 120. );
//    _hzPos         = tfdir.make<TH1F>( "hzPos",
//                                       "MuonicXRay Photon z at Production;(mm)",
//                                       bins.nbins(), bins.low(), bins.high() );
    _htime         = tfdir.make<TH1F>( "htime",
                                       "MuonicXRay Photon time at Production;(ns)",
                                       210, -200., 3000. );
    _hxyPos        = tfdir.make<TH2F>( "hxyPos",
                                       "MuonicXRay Photon (x,y) at Production;(mm)",
                                       60,  -120., 120., 60, -120., 120. );
    //_hrzPos        = tfdir.make<TH2F>( "hrzPos",
    //                                   "MuonicXRay Photon (z,r) at Production;(mm)",
    //                                   bins2.nbins(), bins2.low(), bins2.high(), 60, 0., 120. );

  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::StoppedMuonXRayGammaRayGun);
