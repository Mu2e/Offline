// $Id: ParticleGunImpl.cc,v 1.1 2012/02/03 05:08:06 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/03 05:08:06 $

#include "EventGenerator/inc/ParticleGunImpl.hh"

#include <iostream>

// Framework includes
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Mu2e includes
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"

// Root includes
#include "TH1F.h"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

namespace mu2e {

  ParticleGunImpl::ParticleGunImpl(double meanMultiplicity,
				   PDGCode::type pdgId,
				   double pmin,
				   double pmax,
				   const RandomUnitSphereParams& angles,
				   double tmin,
				   double tmax,

				   const CLHEP::Hep3Vector& point,
				   const CLHEP::Hep3Vector& halfLength,
				   
				   const std::string& histoDir,
				   bool verbose)
    : GeneratorBase()

    , _mean(meanMultiplicity)
    , _pdgId(pdgId)

    , _pmin(pmin)
    , _pmax(pmax)

    , _tmin(tmin)
    , _tmax(tmax)

    , _point(point)

    , _halfLength(halfLength)

    , _verbose(verbose)
    , _doHistograms(!histoDir.empty())

    // Random number distributions; getEngine() comes from base class.
    , _randFlat( getEngine() )
    , _randPoissonQ( getEngine(), std::abs(_mean) )
    , _randomUnitSphere( getEngine(), angles)

    // Histogram pointers
    , _hMultiplicity(0)
    , _hMomentum(0)
    , _hCz(0)
    , _hX0(0)
    , _hY0(0)
    , _hZ0(0)
    , _hT0(0)
  {
    ConditionsHandle<ParticleDataTable> pdt("ignored");
    _mass = pdt->particle(_pdgId).ref().mass().value();
    
    _dp  = ( _pmax - _pmin);
    _dt  = ( _tmax - _tmin);

    // Book histograms if enabled.
    if (_doHistograms) {

      art::ServiceHandle<art::TFileService> tfs;

      art::TFileDirectory tfdir = tfs->mkdir(histoDir);

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
  }

  void ParticleGunImpl::generate( GenParticleCollection& genParts ){

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

      genParts.push_back( GenParticle( _pdgId, GenId::particleGun, pos, p4, time));

      if(_verbose) {
	cout << "Generated position: "<< pos << " "
	     << p4 << " "
	     << p4.vect().mag() << " "
	     << time
	     << endl;
      }

      if ( _doHistograms ) {
        _hMomentum->Fill(p);
        _hCz->Fill( p4.vect().cosTheta());
        _hX0->Fill( pos.x() );
        _hY0->Fill( pos.y() );
        _hZ0->Fill( pos.z() );
        _hT0->Fill( time );
      }
    }
  } // generate()

} // namespace mu2e
