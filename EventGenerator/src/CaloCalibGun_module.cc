// Simulate the photons coming from the pipe calibration source
// based on CaloCalibGun orginally written by Bertrand Echenard (2014)
// Current module author: Sophie Middleton (2022)

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

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
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh" 

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandPoissonQ.h"


#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//ROOT Includes
#include <TMath.h>

const double piconst = CLHEP::pi;
using CLHEP::degree;
using namespace std;
using namespace mu2e;

namespace mu2e {

  class CaloCalibGun : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<unsigned> mean{Name("mean"),1};
      fhicl::Atom<double> energy{Name("PhotonEnergy"),6.13};//MeV
      fhicl::Atom<double> cosmin{Name("cosmin"),-1.};
      fhicl::Atom<double> cosmax{Name("cosmax"),1.};
      fhicl::Atom<double> phimin{Name("phimin"),0.};
      fhicl::Atom<double> phimax{Name("phimax"), CLHEP::twopi };
      fhicl::Atom<double> tmin{Name("tmin"),0.};
      fhicl::Atom<double> tmax{Name("tmax"),1694.};
      fhicl::Atom<bool> doHistograms{Name("doHistograms"), true };
    };

    using Parameters= art::EDProducer::Table<Config>;
    explicit CaloCalibGun(const Parameters& conf);

    virtual void produce(art::Event& event) override;

    //----------------------------------------------------------------
  private:

    unsigned _mean;
    double _energy;
    double _cosmin;
    double _cosmax;
    double _phimin;
    double _phimax;
    double _tmin;
    double _tmax;
    
    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandFlat     _randFlat;
    CLHEP::RandPoissonQ _randPoissonQ;
    RandomUnitSphere    _randomUnitSphere;

    const DiskCalorimeter *_cal;
    //int                    _nPipes;
    double                 _pipeRadius;
    std::vector<double>    _pipeTorRadius;
    std::vector<double>    _randomRad;
    CLHEP::Hep3Vector      _zPipeCenter;

  };

  CaloCalibGun::CaloCalibGun(const Parameters& conf)
    : EDProducer{conf}
    , _mean{conf().mean()}
    , _energy{conf().energy()}
    , _cosmin{conf().cosmin()}
    , _cosmax{conf().cosmax()}
    , _phimin{conf().phimin()}
    , _phimax{conf().phimax()}
    , _tmin{conf().tmin()}
    , _tmax{conf().tmax()}
    , _engine{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , _randFlat{_engine}
    , _randPoissonQ{_engine, _mean*1.0}
    , _randomUnitSphere{_engine, _cosmin, _cosmax, 0, CLHEP::twopi}
  {
    produces<mu2e::GenParticleCollection>();
    produces<mu2e::PrimaryParticle>();
    produces <MCTrajectoryCollection>();

    if (_mean < 0) throw cet::exception("RANGE") << "CaloCalibGun.mean must be non-negative "<< '\n';
  }

  //================================================================
  void CaloCalibGun::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    PrimaryParticle primaryParticles;
    MCTrajectoryCollection mctc;

    _cal =  &*GeomHandle<mu2e::DiskCalorimeter>();
    //_nPipes          = _cal->caloInfo().getInt("nPipes");
    _pipeRadius      = _cal->caloInfo().getDouble("pipeRadius");
    _pipeTorRadius   = _cal->caloInfo().getVDouble("pipeTorRadius");
    _randomRad       = _cal->caloInfo().getVDouble("pipeTorRadius");
    _zPipeCenter     = _cal->disk(1).geomInfo().origin()-CLHEP::Hep3Vector(0,0,_cal->disk(1).geomInfo().size().z()/2.0-_pipeRadius);//disk =1
    // we normalize to the volume of the pipe (proportional to 2*pi*R if they have all the same radius) to draw a random number from which to generate the photons TODO - is this used?
    double sumR(0);
    std::for_each(_randomRad.begin(), _randomRad.end(), [&](double& d) {sumR+=d; d = sumR; });
    std::for_each(_randomRad.begin(), _randomRad.end(), [&](double& d) {d /= sumR;});
    
    unsigned int nGen = _mean;
    //Define the parameters of the pipes:
    
    // angle of large torus in degrees
    std::vector<double> phi_lbd = _cal->caloInfo().getVDouble("LargeTorPhi");//{161.34, 149.50, 139.50, 132.07, 125.39}; 
    // angle of small torus in degrees
    std::vector<double> phi_sbd = _cal->caloInfo().getVDouble("smallTorPhi");//{84.63, 85.28, 85.79, 86.20, 86.53}; 
    // angle of the end point
    std::vector<double> phi_end = _cal->caloInfo().getVDouble("straitEndPhi");//{3.96, 10.53, 15.80, 20.16, 23.84};
    // center position y of the small torus
    std::vector<double> ysmall = _cal->caloInfo().getVDouble("yposition");//{432.2, 480.5, 524.3, 564.7, 602.5}; 
    // radius of small torus
    double radSmTor = _cal->caloInfo().getDouble("radSmTor");//41.0;
    // first center position x of the small torus
    double xsmall = _cal->caloInfo().getDouble("radSmTor");//71.0;
    // distance of the small torus center
    double xdistance = _cal->caloInfo().getDouble("xdistance");//60.0;
    // inner radius of the manifold
    double rInnerManifold = _cal->caloInfo().getDouble("rInnerManifold");//681.6; // 713.35 mm - 1.25 in (31.75 mm)
    vector<float> sign{-1.0, 1.0};

    for (unsigned int ig = 0; ig < nGen; ++ig) 
    {
      // only in the front of 2nd disk.
      //int nDisk        = 1;
      //_zPipeCenter     = _cal->disk(nDisk).geomInfo().origin()-CLHEP::Hep3Vector(0,0,_cal->disk(nDisk).geomInfo().size().z()/2.0-_pipeRadius);

      double xpipe, ypipe, zpipe;
      //Pick position
      int xsn = round(_randFlat.fire());
      int ysn = round(_randFlat.fire());

      double theta = _randFlat.fire() * 2.0 * piconst;
      double pipeR = _pipeRadius * _randFlat.fire();
      zpipe = pipeR*sin(theta);

      double rtest = _randFlat.fire();
      int idx = int( std::lower_bound(_randomRad.begin(), _randomRad.end(), rtest) - _randomRad.begin());
      double radLgTor = _pipeTorRadius[idx];

      // The phi range from 0 to half phi_lbd for the large torus
      double phiLgTor = _randFlat.fire() * phi_lbd[idx]  * CLHEP::degree / 2.;
      // x, y, z position of the large torus
      double xLgTor = sign[xsn]*(radLgTor + pipeR*cos(theta))*cos(phiLgTor);
      double yLgTor = sign[ysn]*(radLgTor + pipeR*cos(theta))*sin(phiLgTor);
      // circulus (center) of the large torus
      double circLgTor = radLgTor * phi_lbd[idx]  * CLHEP::degree / 2.;

      // The phi range for the small torus
      double phiSmTor = CLHEP::degree * (_randFlat.fire() * phi_sbd[idx] + 180. + phi_lbd[idx]/2. - phi_sbd[idx]);
      // x, y, z position of the small torus
      double xSmTor = sign[xsn]*((radSmTor + pipeR*cos(theta))*cos(phiSmTor) + xsmall + xdistance * idx);
      double ySmTor = sign[ysn]*((radSmTor + pipeR*cos(theta))*sin(phiSmTor) + ysmall[idx]);
      // circulus (center) of the small torus
      double circSmTor = CLHEP::degree * radSmTor * phi_sbd[idx];

      // strait pipe
      //double xmanifold = rInnerManifold * cos(piconst * (90 - phi_end[idx])/180.);
      double ymanifold = rInnerManifold * sin(CLHEP::degree * (90 - phi_end[idx]));
      double xstart = xsmall + xdistance * idx - radSmTor * cos(CLHEP::degree * phi_end[idx]);
      double ystart = ysmall[idx] + radSmTor * sin(CLHEP::degree * phi_end[idx]);
      // height of the strait pipe
      double hPipe = (ymanifold - ystart) / sin(CLHEP::degree * (90 - phi_end[idx]));
      // a cylinder along y-axis
      double y_center = _randFlat.fire() * hPipe;
      double xPipe = pipeR * cos(theta);
      double xStrait = sign[xsn] * (xPipe * cos(-CLHEP::degree * phi_end[idx]) - y_center * sin(-CLHEP::degree * phi_end[idx] ) + xstart);
      double yStrait = sign[ysn] * (xPipe * sin(-CLHEP::degree * phi_end[idx]) + y_center * cos(-CLHEP::degree * phi_end[idx]) + ystart);
      double lenStrait = hPipe;

      double sample = _randFlat.fire();
      if(sample <= circLgTor / (circLgTor + circSmTor + lenStrait))
      {
        xpipe = xLgTor;
        ypipe = yLgTor;
      }
      else if(sample > circLgTor / (circLgTor + circSmTor + lenStrait) && sample <= (circLgTor + circSmTor) / (circLgTor + circSmTor + lenStrait))
      {
        xpipe = xSmTor;
        ypipe = ySmTor;
      }
      else
      {
        xpipe = xStrait;
        ypipe = yStrait;
      }
      CLHEP::Hep3Vector pos(xpipe, ypipe, zpipe);
      // shift the pipe to the front of the calorimeter disk 0
      pos +=_zPipeCenter;

      //pick time
      double time = _tmin + _randFlat.fire() * ( _tmax - _tmin );

      //Pick energy and momentum vector
      double energy =_energy;
 
      CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_energy);

      //Set Four-momentum
      CLHEP::HepLorentzVector mom(p3.x(), p3.y(),p3.z(),energy );

      // Add the particle to  the list.
      output->emplace_back(PDGCode::gamma, GenId::CaloCalib, pos, mom, time);
      event.put(std::move(output));
    } 
    event.put(std::make_unique<PrimaryParticle>(primaryParticles));
    event.put(std::make_unique<MCTrajectoryCollection>(mctc));
  }
      
  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloCalibGun);
