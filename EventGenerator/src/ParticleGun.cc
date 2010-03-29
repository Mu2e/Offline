//
// Shoots a single particle gun and puts its output
// into a generated event.
//
// $Id: ParticleGun.cc,v 1.3 2010/03/29 16:20:57 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/03/29 16:20:57 $
//
// Original author Rob Kutschke
// 

#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

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

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::twopi;


namespace mu2e {

  ParticleGun::ParticleGun( edm::Run const& , const SimpleConfig& config ):
    GeneratorBase(),
    _n(1),
    _point(),
    _halfLength(),
    _doHistograms(false),
    _hMultiplicity(0),
    _hMomentum(0){

    // Finish setting up default values.
    _point.setZ( defaultZ0() );

    // Conversion energy for Al.  Should come from conditions.
    static const double pEndPoint = 104.96;

    // Process the run time configuration.
    _n = config.getInt("particleGun.n", _n);

    _pdgId = static_cast<PDGCode::type>(config.getInt("particleGun.id",  PDGCode::e_minus));

    if ( config.hasName("particleGun.point") ){
      _point  = config.getHep3Vector("particleGun.point");
    }

    if ( config.hasName("particleGun.halfLength") ){
      _halfLength = config.getHep3Vector("particleGun.halfLength");
    }

    _pmin   = config.getDouble("particleGun.tmin", pEndPoint );
    _pmax   = config.getDouble("particleGun.tmax", pEndPoint );

    _randomUnitSphere.setczmin(  config.getDouble("particleGun.czmin",  0.5));
    _randomUnitSphere.setczmax(  config.getDouble("particleGun.czmax",  0.7));
    _randomUnitSphere.setphimin( config.getDouble("particleGun.czmin",   0.));
    _randomUnitSphere.setphimax( config.getDouble("particleGun.czmax",  twopi));

    _tmin   = config.getDouble("particleGun.tmin",  0. );
    _tmax   = config.getDouble("particleGun.tmax",  0. );

    _doHistograms = config.getBool("particleGun.doHistograms", _doHistograms);

    // Compute derived information.
    ConditionsHandle<ParticleDataTable> pdt("ignored");
    _mass = pdt->particle(_pdgId).mass().value();

    _dp  = ( _pmax - _pmin);
    _dt  = ( _tmax - _tmin);

    // Book histograms if enabled.
    if ( !_doHistograms ) return;

    edm::Service<edm::TFileService> tfs;

    edm::TFileDirectory tfdir = tfs->mkdir( "ParticleGun" );
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Particle Gun Multiplicity",    10,  0.,  10.);
    _hMomentum     = tfs->make<TH1F>( "hMomentum",     "Particle Gun Momentum (MeV)",  100, 0., 110.);
    _hCz           = tfs->make<TH1F>( "hCz",           "Particle Gun cos(theta)",      100, -1.,  1.);
    _hX0           = tfs->make<TH1F>( "hX0",           "Particle Gun X0",              100,   -20.,    20.);
    _hY0           = tfs->make<TH1F>( "hY0",           "Particle Gun Y0",              100,   -20.,    20.);
    _hZ0           = tfs->make<TH1F>( "hZ0",           "Particle Gun Z0",              100, -2000.,  2000.);
    _hT0           = tfs->make<TH1F>( "hT0",           "Particle Gun Time",            100,     0.,  2000.);

  }

  ParticleGun::~ParticleGun(){
  }

  void ParticleGun::generate( ToyGenParticleCollection& genParts ){

    if ( _doHistograms ){
      _hMultiplicity->Fill(_n);
    }

    for ( int j =0; j<_n; ++j ){

      // Random point inside of a box.
      double x[3];
      for ( int i=0; i<3; ++i ){
        x[i] = _point[i] + RandFlat::shoot(-1.,1.)*_halfLength[i];
      }
      Hep3Vector pos(x[0], x[1], x[2]);

      // Magnitude of momentum.
      double p = _pmin + _dp * RandFlat::shoot();
      
      // 3 Momentum.
      Hep3Vector mom = p*_randomUnitSphere.shoot();

      // Energy
      double e = sqrt(mom*mom + _mass*_mass);
    
      // 4 Momentum.
      HepLorentzVector p4( mom, e);

      // Time
      double time = _tmin + _dt*RandFlat::shoot();

      genParts.push_back( ToyGenParticle( _pdgId, GenId::particleGun, pos, p4, time));

      if ( _doHistograms ) {
        _hMomentum->Fill(p);
        _hCz->Fill(mom.cosTheta());
        _hX0->Fill( pos.x() );
        _hY0->Fill( pos.y() );
        _hZ0->Fill( pos.z() );
        _hT0->Fill( time );
      }

    }

  } // end of generate


  // Compute a default value of the z0 of the particle gun.
  // In the tracker coordinate system, in mm.
  double ParticleGun::defaultZ0(){

    // Default z position is at the entrance to a nominal tracker, in mm.
    double z0 = -1400.;

    // Override this if an actual tracker is present.

    edm::Service<GeometryService> geom;

    if ( geom->hasElement<LTracker>() ){
      GeomHandle<LTracker> ltracker;
      z0 = -ltracker->zHalfLength();
    } else if ( geom->hasElement<ITracker>() ) {
      GeomHandle<ITracker> itracker;
      z0 = -itracker->zHalfLength();
    }

    // Add a TTracker branch when ready.

    // Tracker position in the Detecctor coordinate system.
    static double trackerOffset = -1800.;

    return z0+trackerOffset;
  }

}
