//
//
// Simulate the photons coming from the pipe calibration source
// //
//
// Original author Bertrand Echenard
//
//

// C++ includes.
#include <iostream>
#include <algorithm>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/CaloCalibGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

using namespace std;

namespace mu2e {

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

    // Histogram pointers
    _hE(0),_hT(0),_hcos(0),_hphi(0),_hrad(0),_hz(0),_hxy(0)
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
    _zPipeCenter     = _cal->disk(0).geomInfo().origin()-CLHEP::Hep3Vector(0,0,_cal->disk(0).geomInfo().size().z()/2.0-_pipeRadius);




    // we normalize to the volume of the pipe (proportional to 2*pi*R if they have all the same radius) to draw a
    // random number from which to generate the photons
    double sumR(0);
    std::for_each(_randomRad.begin(), _randomRad.end(), [&](double& d) {sumR+=d; d = sumR; });
    std::for_each(_randomRad.begin(), _randomRad.end(), [&](double& d) {d /= sumR;});

    // Book histograms.
    if ( _doHistograms )
    {
       art::ServiceHandle<art::TFileService> tfs;
       art::TFileDirectory tfdir  = tfs->mkdir( "CaloPhotonCalibGun" );
       _hE   = tfdir.make<TH1D>( "hE",   "Photon Energy"   ,  100,      0, 20);
       _hT   = tfdir.make<TH1D>( "hT",   "Photon Time"     ,  100,      0, 2000);
       _hcos = tfdir.make<TH1D>( "hcos", "Photon cos theta",  100,     -1, 1);
       _hphi = tfdir.make<TH1D>( "hphi", "Photon phi"      ,  100,      0, 6.3);
       _hrad = tfdir.make<TH1D>( "hrad", "Pos radius"      ,  100,      0, 700);
       _hz   = tfdir.make<TH1D>( "hz",   "Pos z "          ,  100,  10000, 14000);
       _hxy  = tfdir.make<TH2D>( "hxy",  "Pos xy"          ,  100,   -700, 700, 100, -700 ,700);
    }


    if (_mean < 0) throw cet::exception("RANGE") << "CaloCalibGun.mean must be non-negative "<< '\n';

  }


// to do on that one
// get the correct z position of the pipes
// check what these detector positions are



  CaloCalibGun::~CaloCalibGun(){}



  void CaloCalibGun::generate( GenParticleCollection& genParts ){


      //int nGen = _randPoissonQ.fire();
      int nGen = _mean;
      for (int ig=0; ig<nGen; ++ig) {

        //Pick position
        double rtest = _randFlat.fire();
        int idx = int( std::lower_bound(_randomRad.begin(), _randomRad.end(), rtest) - _randomRad.begin());
        double rad = _pipeTorRadius[idx];
        double phi = _randFlat.fire()*(_phimax-_phimin)+_phimin;

        CLHEP::Hep3Vector pos(rad*cos(phi),rad*sin(phi),0);
        pos +=_zPipeCenter;



        //pick time
        double time = _tmin + _randFlat.fire() * ( _tmax - _tmin );

        //Pick energy and momentum vector
        //double e = _elow + _flatmomentum.fire() * ( _ehi - _elow );
        double energy =_energy;
        CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_energy);
        while (p3.cosTheta()<0) p3 = _randomUnitSphere.fire(_energy);

        //Set Four-momentum
        CLHEP::HepLorentzVector mom(0,0,0,0);
        mom.setPx( p3.x() );
        mom.setPy( p3.y() );
        mom.setPz( p3.z() );
        mom.setE( energy );

        // Add the particle to  the list.
        genParts.push_back( GenParticle(PDGCode::gamma, GenId::CaloCalib, pos, mom, time));

        // Fill histograms.
        if ( _doHistograms) {
          _hE    ->Fill(energy);
          _hT    ->Fill(time);
          _hcos  ->Fill(p3.cosTheta());
          _hphi  ->Fill(phi);
          _hrad  ->Fill(rad);
          _hz    ->Fill(pos.z());
          _hxy   ->Fill(rad*cos(phi),rad*sin(phi));
        }
      }


  }

}
