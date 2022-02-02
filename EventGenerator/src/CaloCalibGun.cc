//
// Simulate the photons coming from the pipe calibration source
// //
// $Id: CaloCalibGun.cc,v 1.16 2014/01/27 22:20:17 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/01/27 22:20:17 $
//
// Original author Bertrand Echenard
// Edited by: De Xu Lin, Sophie Middleton

// C++ includes.
#include <iostream>
#include <algorithm>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/ConditionsService/inc/AcceleratorParams.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataTable.hh"
#include "Offline/EventGenerator/inc/CaloCalibGun.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TMath.h"
using namespace std;
using namespace TMath;
const double piconst = Pi();

namespace mu2e
{
  CaloCalibGun::CaloCalibGun(CLHEP::HepRandomEngine& engine, art::Run& run, const SimpleConfig& config):
  // Configurable parameters
  _mean(config.getDouble("caloCalibGun.mean",1.)),
  _energy(config.getDouble("caloCalibGun.energy",6.13)),
  _cosmin(config.getDouble("caloCalibGun.cosmin",  -1.)),
  _cosmax(config.getDouble("caloCalibGun.cosmax",  1.)),
  _phimin(config.getDouble("caloCalibGun.phimin", 0. )),
  _phimax(config.getDouble("caloCalibGun.phimax", CLHEP::twopi )),
  _randFlat{engine},
  _randPoissonQ{engine, std::abs(_mean)},
  _randomUnitSphere{engine, _cosmin, _cosmax, 0, CLHEP::twopi},
  _detSys(),
  _doHistograms(config.getBool("caloCalibGun.doHistograms",true)),
  _Ntupe(0)
  {
    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    GlobalConstantsHandle<ParticleDataTable> pdt;

    // Default values for the start and end of the live window.
    _tmin = 0.;
    _tmax = accPar->deBuncherPeriod;
    _tmin = config.getDouble("caloCalibGun.tmin",  _tmin );
    _tmax = config.getDouble("caloCalibGun.tmax",  _tmax );

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
      _Ntupe  = tfs->make<TTree>("caliGun", "caliGun");

      _Ntupe -> Branch("genErg", &_genErg, "genErg/F");
      _Ntupe -> Branch("genTime", &_genTime, "genTime/F");
      _Ntupe -> Branch("genCos", &_genCos, "genCos/F");
      //_Ntupe -> Branch("genPhi", &_genPhi, "genPhi/F");
      //_Ntupe -> Branch("genRad", &_genRad, "genRad/F");
      _Ntupe -> Branch("genPosX", &_genPosX, "genPosX/F");
      _Ntupe -> Branch("genPosY", &_genPosY, "genPosY/F");
      _Ntupe -> Branch("genPosZ", &_genPosZ, "genPosZ/F");
    }
    if (_mean < 0) throw cet::exception("RANGE") << "CaloCalibGun.mean must be non-negative "<< '\n';
  }


  // to do on that one
  // get the correct z position of the pipes
  // check what these detector positions are
  CaloCalibGun::~CaloCalibGun(){}
  void CaloCalibGun::generate(GenParticleCollection& genParts)
  {
    //int nGen = _randPoissonQ.fire();
    int nGen = _mean;
    // define the parameters of the pipes
    // angle of large torus in degrees
    double phi_lbd[5] = {161.34, 149.50, 139.50, 132.07, 125.39};
    // angle of small torus in degrees
    double phi_sbd[5] = {84.63, 85.28, 85.79, 86.20, 86.53};
    // angle of the end point
    double phi_end[5] = {3.96, 10.53, 15.80, 20.16, 23.84};
    // center position y of the small torus
    double ysmall[5] = {432.2, 480.5, 524.3, 564.7, 602.5};
    // radius of small torus
    double radSmTor = 41.0;
    // first center position x of the small torus
    double xsmall = 71.0;
    // distance of the small torus center
    double xdistance = 60.0;
    // inner radius of the manifold
    double rInnerManifold = 681.6; // 713.35 mm - 1.25 in (31.75 mm)
    vector<float> sign{-1.0, 1.0};

    for (int ig=0; ig<nGen; ++ig) 
    {
      // pick up the disk, 0 or 1 
      //int nDisk        = rint(_randFlat.fire());
      // only in the front of 2nd disk.
      int nDisk        = 1;
      _zPipeCenter     = _cal->disk(nDisk).geomInfo().origin()-CLHEP::Hep3Vector(0,0,_cal->disk(nDisk).geomInfo().size().z()/2.0-_pipeRadius);
      //cout << "The number of disk: " << nDisk << "; z position: " << _zPipeCenter << endl;
 
      double xpipe, ypipe, zpipe;
      //Pick position
      int xsn = rint(_randFlat.fire());
      int ysn = rint(_randFlat.fire());

      double theta = _randFlat.fire() * 2.0 * piconst;
      double pipeR = _pipeRadius * _randFlat.fire();
      zpipe = pipeR*sin(theta);

      double rtest = _randFlat.fire();
      int idx = int( std::lower_bound(_randomRad.begin(), _randomRad.end(), rtest) - _randomRad.begin());
      double radLgTor = _pipeTorRadius[idx];
      // modify to the full circle, Aug. 21, 2018
      // The phi range from 0 to half phi_lbd for the large torus
      double phiLgTor = _randFlat.fire() * phi_lbd[idx]  * piconst / 2. / 180.;
      // x, y, z position of the large torus
      double xLgTor = sign[xsn]*(radLgTor + pipeR*cos(theta))*cos(phiLgTor);
      double yLgTor = sign[ysn]*(radLgTor + pipeR*cos(theta))*sin(phiLgTor);
      // circulus (center) of the large torus
      double circLgTor = radLgTor * phi_lbd[idx]  * piconst / 2. / 180.;

      // The phi range for the small torus
      double phiSmTor = piconst * (_randFlat.fire() * phi_sbd[idx] + 180. + phi_lbd[idx]/2. - phi_sbd[idx])/180.;
      // x, y, z position of the small torus
      double xSmTor = sign[xsn]*((radSmTor + pipeR*cos(theta))*cos(phiSmTor) + xsmall + xdistance * idx);
      double ySmTor = sign[ysn]*((radSmTor + pipeR*cos(theta))*sin(phiSmTor) + ysmall[idx]);
      // circulus (center) of the small torus
      double circSmTor = radSmTor * piconst * phi_sbd[idx] / 180.;

      // strait pipe
      //double xmanifold = rInnerManifold * cos(piconst * (90 - phi_end[idx])/180.);
      double ymanifold = rInnerManifold * sin(piconst * (90 - phi_end[idx])/180.);
      double xstart = xsmall + xdistance * idx - radSmTor * cos(piconst * phi_end[idx]/180.);
      double ystart = ysmall[idx] + radSmTor * sin(piconst * phi_end[idx]/180.);
      // height of the strait pipe
      double hPipe = (ymanifold - ystart) / sin(piconst * (90 - phi_end[idx]) / 180.);
      // a cylinder along y-axis
      double y_center = _randFlat.fire() * hPipe;
      double xPipe = pipeR * cos(theta);
      double xStrait = sign[xsn] * (xPipe * cos(-piconst * phi_end[idx] / 180.) - y_center * sin(-piconst * phi_end[idx] / 180.) + xstart);
      double yStrait = sign[ysn] * (xPipe * sin(-piconst * phi_end[idx] / 180.) + y_center * cos(-piconst * phi_end[idx] / 180.) + ystart);
      double lenStrait = hPipe;

      double sample = _randFlat.fire();
      if(sample <= circLgTor / (circLgTor + circSmTor + lenStrait))
//      if(sample <= circLgTor / (circLgTor + circSmTor))
      {
        xpipe = xLgTor;
        ypipe = yLgTor;
      }
      else if(sample > circLgTor / (circLgTor + circSmTor + lenStrait) && sample <= (circLgTor + circSmTor) / (circLgTor + circSmTor + lenStrait))
//      else
      {
        xpipe = xSmTor;
        ypipe = ySmTor;
      }
      else
      {
        xpipe = xStrait;
        ypipe = yStrait;
      }

      //cout << "===============================================" << endl;
      //cout << "pipe center position: " << _zPipeCenter.x() << "; " << _zPipeCenter.y() << "; " << _zPipeCenter.z() << endl;
      //cout << "X position: " << rad*cos(phi) << endl;
      CLHEP::Hep3Vector pos(xpipe, ypipe, zpipe);
      // shift the pipe to the front of the calorimeter disk 0
      pos +=_zPipeCenter;
      //cout << "Positions: x = " << pos.x() << "; y = " << pos.y() << "; z = " << pos.z() << endl;

      //pick time
      double time = _tmin + _randFlat.fire() * ( _tmax - _tmin );

      //Pick energy and momentum vector
      //double e = _elow + _flatmomentum.fire() * ( _ehi - _elow );
      double energy =_energy;
      CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_energy);
      //while (p3.cosTheta()<0) p3 = _randomUnitSphere.fire(_energy);

      //Set Four-momentum
      CLHEP::HepLorentzVector mom(0,0,0,0);
      mom.setPx( p3.x() );
      mom.setPy( p3.y() );
      mom.setPz( p3.z() );
      mom.setE( energy );

      // Add the particle to  the list.
      genParts.push_back( GenParticle(PDGCode::gamma, GenId::CaloCalib, pos, mom, time));

      // Fill histograms.
      if ( _doHistograms)
      {
        _genErg = energy;
        _genTime = time;
        _genCos = p3.cosTheta();
        //_genPhi = phi;
        //_genRad = rad;
        _genPosX = pos.x();
        _genPosY = pos.y();
        _genPosZ = pos.z();
      }
      if(_doHistograms) _Ntupe -> Fill();
    }

    
  }
}
