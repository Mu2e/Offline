//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: ParticleGun.cc,v 1.11 2011/05/18 14:21:44 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/18 14:21:44 $
//
// Original author Rob Kutschke
//

#include <iostream>

// Framework includes
#include "art/Framework/Core/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/ParticleGun.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// Root includes
#include "TH1F.h"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

//using CLHEP::Hep3Vector;
//using CLHEP::HepLorentzVector;
//using CLHEP::RandFlat;
//using CLHEP::twopi;


namespace mu2e {

  // Conversion energy for Al.  Should come from conditions.
  static const double pEndPoint = 104.96;

  ParticleGun::ParticleGun( art::Run const& , const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // From run time configuration file.
    _mean(config.getDouble("particleGun.mean",-1.)),
    _pdgId(static_cast<PDGCode::type>(config.getInt("particleGun.id",  PDGCode::mu_minus))),
    _czmin( config.getDouble("particleGun.czmin",  0.5)),
    _czmax( config.getDouble("particleGun.czmax",  0.7)),
    _phimin(config.getDouble("particleGun.phimin", 0. )),
    _phimax(config.getDouble("particleGun.phimax", CLHEP::twopi)),
    _pmin(pEndPoint),
    _pmax(pEndPoint),
    _tmin(0),
    _tmax(0),
    _point(),
    _halfLength(),
    _doHistograms(config.getBool("particleGun.doHistograms", false)),

    // Random number distributions; getEngine() comes from base class.
    _randFlat( getEngine() ),
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randomUnitSphere( getEngine(), _czmin, _czmax, _phimin, _phimax),

    // Histogram pointers
    _hMultiplicity(0),
    _hMomentum(0),
    _hCz(0),
    _hX0(0),
    _hY0(0),
    _hZ0(0),
    _hT0(0){

    if ( config.hasName("particleGun.point") ){
      _point  = config.getHep3Vector("particleGun.point");
    }

    if ( config.hasName("particleGun.halfLength") ){
      _halfLength = config.getHep3Vector("particleGun.halfLength");
    }

    _pmin   = config.getDouble("particleGun.pmin", pEndPoint );
    _pmax   = config.getDouble("particleGun.pmax", pEndPoint );

    _tmin   = config.getDouble("particleGun.tmin",  0. );
    _tmax   = config.getDouble("particleGun.tmax",  0. );

    // end processing run time configuration.
    ConditionsHandle<ParticleDataTable> pdt("ignored");
    _mass = pdt->particle(_pdgId).ref().mass().value();

    _dp  = ( _pmax - _pmin);
    _dt  = ( _tmax - _tmin);

    // Book histograms if enabled.
    if ( !_doHistograms ) return;

    art::ServiceHandle<art::TFileService> tfs;

    art::TFileDirectory tfdir = tfs->mkdir( "ParticleGun" );
    _hMultiplicity = tfdir.make<TH1F>( "hMultiplicity", "Particle Gun Multiplicity",    20,  0.,  20.);

    // Pick range of histogram.  Nice round numbers for the usual case.
    if ( _pmax > 100. && _pmax < 107. ){
      _hMomentum     = tfdir.make<TH1F>( "hMomentum",     "Particle Gun Momentum (MeV)",  100, 0., 110.);
    } else {
      double pHigh = double(static_cast<int>(_pmax)+1);
      _hMomentum     = tfdir.make<TH1F>( "hMomentum",     "Particle Gun Momentum (MeV)",  100, 0., pHigh );
    }

    _hCz           = tfdir.make<TH1F>( "hCz",           "Particle Gun cos(theta)",      100, -1.,  1.);

    double xlen = (_halfLength.x() > 1. ) ? _halfLength.x() : 1.;
    double ylen = (_halfLength.y() > 1. ) ? _halfLength.y() : 1.;
    double zlen = (_halfLength.z() > 1. ) ? _halfLength.z() : 1.;
    double tlen = ((_tmax-_tmin)   > 1. ) ? (_tmax-_tmin)   : 1.;
    double xl = _point.x() - xlen;
    double xh = _point.x() + xlen;
    double yl = _point.y() - ylen;
    double yh = _point.y() + ylen;
    double zl = _point.z() - zlen;
    double zh = _point.z() + zlen;
    _hX0           = tfdir.make<TH1F>( "hX0", "Particle Gun X0",              100,   xl,    xh);
    _hY0           = tfdir.make<TH1F>( "hY0", "Particle Gun Y0",              100,   yl,    yh);
    _hZ0           = tfdir.make<TH1F>( "hZ0", "Particle Gun Z0",              100,   zl,    zh);
    _hT0           = tfdir.make<TH1F>( "hT0", "Particle Gun Time",            100,  _tmin, _tmin + tlen);


  }

  ParticleGun::~ParticleGun(){
  }

  void ParticleGun::generate( ToyGenParticleCollection& genParts ){

    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ){
      _hMultiplicity->Fill(n);
    }

    for ( int j =0; j<n; ++j ){

      // Random point inside of a box.
      double x[3];
      for ( int i=0; i<3; ++i ){
        x[i] = _point[i] + _randFlat.fire(-1.,1.)*_halfLength[i];
      }
      CLHEP::Hep3Vector pos(x[0], x[1], x[2]);

      // Magnitude of momentum and energy.
      double p = _pmin + _dp * _randFlat.fire();
      double e = sqrt( p*p + _mass*_mass);

      // 4 Momentum.
      CLHEP::HepLorentzVector p4( _randomUnitSphere.fire(p), e);

      // Time
      double time = _tmin + _dt*_randFlat.fire();

      genParts.push_back( ToyGenParticle( _pdgId, GenId::particleGun, pos, p4, time));

      cout << "Generated position: "
           << pos << " "
           << p4 << " "
           << p4.vect().mag() << " "
           << time
           << endl;

      if ( _doHistograms ) {
        _hMomentum->Fill(p);
        _hCz->Fill( p4.vect().cosTheta());
        _hX0->Fill( pos.x() );
        _hY0->Fill( pos.y() );
        _hZ0->Fill( pos.z() );
        _hT0->Fill( time );
      }

    }

  } // end of generate

}
