// $Id: ParticleGunImpl.cc,v 1.6 2013/03/22 22:04:54 gandr Exp $
// $Author: gandr $
// $Date: 2013/03/22 22:04:54 $
// Original author Robert Kutschke
// Modified by mjlee. See docdb-2049

#include "EventGenerator/inc/ParticleGunImpl.hh"

#include <iostream>

// Framework includes
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"

// Root includes
#include "TH1F.h"
#include "TNtuple.h"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include <boost/algorithm/string.hpp>

using namespace std;

namespace mu2e {

  ParticleGunImpl::ParticleGunImpl(
      CLHEP::HepRandomEngine& engine,
      double meanMultiplicity,
      PDGCode::type pdgId,
      double pmin,
      double pmax,
      const std::string& momentumModeString,

      const RandomUnitSphereParams& angles,
      double tmin,
      double tmax,

      const CLHEP::Hep3Vector& point,
      const CLHEP::Hep3Vector& halfLength,
      const std::string& sourceShapeString,

      int  iterationLimit,
      bool throwOnIterationLimit,

      bool useDetectorCoordinateSystem,

      const std::string& histoDir,
      bool doNtuples,

      int verbosityLevel)
    : _mean(meanMultiplicity)
    , _pdgId(pdgId)

    , _pmin(pmin)
    , _pmax(pmax)

    , _tmin(tmin)
    , _tmax(tmax)

    , _point(point)

    , _halfLength(halfLength)

    , _iterationLimit(iterationLimit)
    , _throwOnIterationLimit(throwOnIterationLimit)

    , _doNtuples(doNtuples)

    , _verbosityLevel(verbosityLevel)
    , _doHistograms(!histoDir.empty())
    , _histoDir(histoDir)
    , _useDetectorCoordinateSystem(useDetectorCoordinateSystem)

    , _randFlat( engine )
    , _randPoissonQ( engine, std::abs(_mean) )
    , _randGaussQ( engine )
    , _randExponential ( engine )
    , _randomUnitSphere( engine, angles)

    // Histogram pointers
    , _hMultiplicity(0)
    , _hMomentum(0)
    , _hCz(0)
    , _hX0(0)
    , _hY0(0)
    , _hZ0(0)
    , _hT0(0)

    , _origin(0)
    , _originSetFlag (false)
  {
    initialize(momentumModeString, sourceShapeString);
  }


  // Old style constructor. Kept for backward compatibility.
  ParticleGunImpl::ParticleGunImpl(
      CLHEP::HepRandomEngine& engine,
      double meanMultiplicity,
      PDGCode::type pdgId,
      double pmin,
      double pmax,

      const RandomUnitSphereParams& angles,
      double tmin,
      double tmax,

      const CLHEP::Hep3Vector& point,
      const CLHEP::Hep3Vector& halfLength,

      const std::string& histoDir,

      int verbosityLevel)
    : _mean(meanMultiplicity)
    , _pdgId(pdgId)

    , _pmin(pmin)
    , _pmax(pmax)

    , _tmin(tmin)
    , _tmax(tmax)

    , _point(point)

    , _halfLength(halfLength)

    , _iterationLimit(100)
    , _throwOnIterationLimit(false)
    , _doNtuples(false)

    , _verbosityLevel(verbosityLevel)
    , _doHistograms(!histoDir.empty())
    , _histoDir (histoDir)
    , _useDetectorCoordinateSystem (false)

    // Random number distributions
    , _randFlat( engine )
    , _randPoissonQ( engine, std::abs(_mean) )
    , _randGaussQ( engine )
    , _randExponential ( engine )
    , _randomUnitSphere( engine, angles)

    // Histogram pointers
    , _hMultiplicity(0)
    , _hMomentum(0)
    , _hCz(0)
    , _hX0(0)
    , _hY0(0)
    , _hZ0(0)
    , _hT0(0)

    , _origin(0)
    , _originSetFlag (false)
  {
    initialize("flat", "box");
  }



  void ParticleGunImpl::initialize (
      const std::string& momentumModeString,
      const std::string& sourceShapeString)
  {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    _mass = pdt->particle(_pdgId).ref().mass().value();

    _dp  = ( _pmax - _pmin);
    _dt  = ( _tmax - _tmin);

    // momentumMode
    std::string momentumModeLower = boost::to_lower_copy(momentumModeString);
    if (momentumModeLower == "flat" ) {
      _momentumMode = flat;
    }
    else if (momentumModeLower == "gaussian") {
      _momentumMode = gaussian;
    }
    else if (momentumModeLower == "generic") {
      _momentumMode = generic;
    }
    else if (momentumModeLower == "exponential") {
      _momentumMode = exponential;
    }
    else if (momentumModeLower == "histogram") {
      _momentumMode = histogram;
    }
    else {
      throw cet::exception("RANGE")
        << "ParticleGun have read invalid momentumMode : "
        << momentumModeString
        << "\n";
    }

    // sourceShape
    std::string sourceShapeLower = boost::to_lower_copy(sourceShapeString);
    if (sourceShapeLower == "box") {
      _sourceShape = box;
    }
    else if (sourceShapeLower == "cylinder_x") {
      _sourceShape = cylinder_x;
    }
    else if (sourceShapeLower == "cylinder_y") {
      _sourceShape = cylinder_y;
    }
    else if (sourceShapeLower == "cylinder_z" || sourceShapeLower == "cylinder") {
      _sourceShape = cylinder_z;
    }
    else if (sourceShapeLower == "sphere") {
      _sourceShape = sphere;
    }
    else {
      throw cet::exception("RANGE")
        << "ParticleGun have read invalid sourceShape : "
        << sourceShapeString
        << "\n";
    }

    // Square of halfLength
    _dx2 = _halfLength.x() * _halfLength.x();
    _dy2 = _halfLength.y() * _halfLength.y();
    _dz2 = _halfLength.z() * _halfLength.z();
    if (    ( _sourceShape == sphere  &&  (_dx2 <=0 || _dy2<=0 || _dz2 <=0) )
         || ( _sourceShape == cylinder_x  &&  (_dy2 <=0 || _dz2<=0 ) )
         || ( _sourceShape == cylinder_y  &&  (_dx2 <=0 || _dz2<=0 ) )
         || ( _sourceShape == cylinder_z  &&  (_dx2 <=0 || _dy2<=0 ) )   ) {

      throw cet::exception("RANGE")
        << "ParticleGun have read strange halfLength values : "
        << _halfLength
        << "\n";
    }


    // Book histograms if enabled.
    if (_doHistograms) {

      art::ServiceHandle<art::TFileService> tfs;

      art::TFileDirectory tfdir = tfs->mkdir(_histoDir);

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

      //Generation Ntuple
      if (_doNtuples) {
        _hGenerationNtuple = tfdir.make<TNtuple>("hGenerationNtuple", "Generation info", "x:y:z:t:px:py:pz:p:pdgId");
      }
    }

    if(_verbosityLevel >= 1) {
      std::cout<<"ParticleGunImpl initialized with mean="<<_mean <<", tmax="<<_tmax<<std::endl;
    }
  }



  void ParticleGunImpl::setParameters (std::vector<double>& momentumParameters) {

/* Here, the parameters for momentum distributions are initialized.
 * pmin and pmax are not initialized here. They should be supplied independently.
 * 1. Flat dist : do not require any parameters. The output will be flat between pmin and pmax.
 * 2. Gaussian dist : par[0] = mean and par[1] = stdDev. After generation, it is checked whether
 *   the result momentum is laid in [pmin, pmax]. Pass-and-fail is applied.
 * 3. Generic dist : More than 2 pars should be supplied. They are treated as a histogram with
 *   exact binning between pmin and pmax. Pass-andfail is applied. No negative or zero parameter allowed.
 * 4. Exponential dist : par[0] is mean of RandExponential, i.e. inverse of decay constant.
 *   pmin is added to the resulting distribution, and pmax boundary check is applied.
 *   Therefore, resulting dist. is f(x) = A exp(- par[0] * (x - pmin)).
 *   par[0] should be >0 .
 * */

    _momentumParameters = momentumParameters;

    if (_momentumMode == gaussian ) {
      if (_momentumParameters.size() !=2) {
        throw cet::exception("RANGE")
          << "ParticleGun require 2 parameters and have read "
          << _momentumParameters.size()
          << " parameters for momentumMode = "
          << _momentumMode
          << "\n";
      }
    }

    else if (_momentumMode == generic) {
      if (_momentumParameters.size() <=1) {
        throw cet::exception("RANGE")
          << "ParticleGun require more than 2 parameters and have read "
          << _momentumParameters.size()
          << " parameters for momentumMode = "
          << _momentumMode
          << "\n";
      }

      double sum = 0;
      for (unsigned int i = 0 ; i <_momentumParameters.size() ; ++i) {
        if (_momentumParameters[i] <0) {
          throw cet::exception("RANGE")
            << "ParticleGun have read negative parameter for momentumMode = "
            << _momentumMode << "\n";
        }
        sum += _momentumParameters[i];
      }

      if (sum <=0) {
        throw cet::exception("RANGE")
          << "ParticleGun have read <=0 sum of momentumParameters for momentumMode = "
          << _momentumMode << "\n";
      }

      for (unsigned int i = 0 ; i <_momentumParameters.size() ; ++i) {
        _momentumParameters[i] /= sum;
      }
    }

    else if (_momentumMode == exponential) {
      if (_momentumParameters.size() < 1) {
        throw cet::exception("RANGE")
          << "ParticleGun require 1 parameters and have read "
          << _momentumParameters.size()
          << " parameters for momentumMode = "
          << _momentumMode
          << "\n";
      }
      if (_momentumParameters[0] <=0) {
        throw cet::exception("RANGE")
          << "ParticleGun have read <=0 mean for momentumMode =  "
          << _momentumMode
          << "\n";
      }
    }
    else if (_momentumMode == histogram) {
      if (_momentumParameters.size() < 3) {
        throw cet::exception("RANGE")
          << "ParticleGun require 3 parameters and have read "
          << _momentumParameters.size()
          << " parameters for momentumMode = "
          << _momentumMode
          << "\n";
      }
      // interpret the 1st 2 parameters as the range
      histogramRange[0] = _momentumParameters[0];
      histogramRange[1] = _momentumParameters[1];
      // the rest are bin contents: normalize them and integrate
      double bint(0.0);
      for(size_t ibin=2;ibin<_momentumParameters.size();++ibin){
        bint += _momentumParameters[ibin];
      }
      double psum(0.0);
      for(size_t ibin=2;ibin<_momentumParameters.size();++ibin){
        histogramBins.push_back(psum + _momentumParameters[ibin]/bint);
        psum = histogramBins.back();
      }
    }

    else {}

  }


  void ParticleGunImpl::generate( GenParticleCollection& genParts ){
    if(_verbosityLevel >= 2) {
      std::cout<<"ParticleGunImpl::generate() for mean="<<_mean <<", tmax="<<_tmax<<std::endl;
    }

    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ){
      _hMultiplicity->Fill(n);
    }
    int iteration (0);

    for ( int j =0; j<n; ++j ){

      // Random point inside of a box.
      double x[3];
      bool passed;

      iteration = 0;
      do {
        passed = true;
        for ( int i=0; i<3; ++i ){
          x[i] = _randFlat.fire(-1.,1.)*_halfLength[i];
        }
        if (_sourceShape == box) break;
        else if (_sourceShape == sphere) {
          if (x[0]*x[0]/_dx2 + x[1]*x[1]/_dy2 + x[2]*x[2]/_dz2 > 1.) passed = false;
        }
        else if (_sourceShape == cylinder_z) {
          if (x[0]*x[0]/_dx2 + x[1]*x[1]/_dy2 >1.) passed = false;
        }
        else if (_sourceShape == cylinder_y) {
          if (x[0]*x[0]/_dx2 + x[2]*x[2]/_dz2 >1.) passed = false;
        }
        else if (_sourceShape == cylinder_x) {
          if (x[1]*x[1]/_dy2 + x[2]*x[2]/_dz2 >1.) passed = false;
        }
        else {}
        ++iteration;
      } while (!passed &&
                   ( (_throwOnIterationLimit && iteration < _iterationLimit) || (!_throwOnIterationLimit ) )
              );

      if (iteration >= _iterationLimit) {
        if (_throwOnIterationLimit) {
          throw cet::exception("RANGE")
            << "ParticleGun have encountered larger number of iteration for vertex generation.\n";
        }
        else {
          mf::LogWarning ("CONTROL")
            << "Warning: ParticleGun have encountered larger number of iteration for vertex generation ("
            << iteration << " / " << _iterationLimit << "\n";
        }

      }

      CLHEP::Hep3Vector pos(x[0]+_point[0], x[1]+_point[1], x[2]+_point[2]);

      // Magnitude of momentum and energy.
      double p(-9999);
      switch (_momentumMode) {

        case flat :
          p = getFlatMomentum();
          break;

        case gaussian :
          p = getGaussianMomentum();
          break;

        case generic :
          p = getGenericMomentum();
          break;

        case exponential :
          p = getExponentialMomentum() ;
          break;

        case histogram :
          p = getHistogramMomentum();
          break;

        default: // Generate flat dist for all unknown momentumMode ?
          p = getFlatMomentum();
          break;
      }

      double e = sqrt( p*p + _mass*_mass);

      // 4 Momentum.
      CLHEP::HepLorentzVector p4( _randomUnitSphere.fire(p), e);

      // Time
      double time = _tmin + _dt*_randFlat.fire();

      // transform to detector coordinate system
      if (_useDetectorCoordinateSystem) {
        if (!_originSetFlag) {
          GeomHandle<DetectorSystem> det;
          _origin = det->toMu2e(CLHEP::Hep3Vector(0.,0.,0.) );
          _originSetFlag = true;
        }
        pos += _origin;
      }

      genParts.push_back( GenParticle( _pdgId, GenId::particleGun, pos, p4, time));

      if(_verbosityLevel >= 3) {
        cout  << "Generated position: "<< pos << " "
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
        if (_doNtuples) {
          _hGenerationNtuple->Fill(pos.x(), pos.y(), pos.z(), time, p4.x(), p4.y(), p4.z(), p4.vect().mag(), _pdgId);
        }
      }
    }
  } // generate()


  double ParticleGunImpl::getFlatMomentum(void) {
    return _pmin + _dp * _randFlat.fire();
  }

  double ParticleGunImpl::getGaussianMomentum(void) {
    int iteration = 0;
    double p(-999999);

    do {
      p = _randGaussQ.fire(_momentumParameters[0], std::abs(_momentumParameters[1]));
      ++iteration;
    } while (  (p < _pmin || p > _pmax) &&
                   ( (_throwOnIterationLimit && iteration < _iterationLimit) || (!_throwOnIterationLimit ) )
            );

    if (iteration >= _iterationLimit) {
      if (_throwOnIterationLimit) {
        throw cet::exception("RANGE")
          << "ParticleGun have encountered larger number of iteration for momentum generation.\n";
      }
      else {
        mf::LogWarning ("CONTROL")
          << "Warning: ParticleGun have encountered larger number of iteration for momentum generation ("
          << iteration << " / " << _iterationLimit << "\n";
      }
    }
    return p;
  }

  double ParticleGunImpl::getGenericMomentum(void) {
    int iteration (0);
    int iindex(0);
    double findex(0);
    double randval(0);
    double idxlim = double(_momentumParameters.size());

    do {
      findex = _randFlat.fire(0., idxlim);
      iindex = int(findex);
      randval = _randFlat.fire();
      ++iteration;
    } while (randval > _momentumParameters[iindex] &&
                   ( (_throwOnIterationLimit && iteration < _iterationLimit) || (!_throwOnIterationLimit ) )
            );

    if (iteration >= _iterationLimit) {
      if (_throwOnIterationLimit) {
        throw cet::exception("RANGE")
          << "ParticleGun have encountered larger number of iteration for momentum generation.\n";
      }
      else {
        mf::LogWarning ("CONTROL")
          << "Warning: ParticleGun have encountered larger number of iteration for momentum generation ("
          << iteration << " / " << _iterationLimit << "\n";
      }
    }
    return _pmin + _dp * findex / idxlim;
  }

  double ParticleGunImpl::getExponentialMomentum(void) {
    int iteration = 0;
    double p(-999999);

    do {
      p = _pmin + _randExponential.fire(_momentumParameters[0]);
      ++iteration;
    } while ((p < _pmin || p > _pmax) &&
                   ( (_throwOnIterationLimit && iteration < _iterationLimit) || (!_throwOnIterationLimit ) )
            );

    if (iteration >= _iterationLimit) {
      if (_throwOnIterationLimit) {
        throw cet::exception("RANGE")
          << "ParticleGun have encountered larger number of iteration for momentum generation.\n";
      }
      else {
        mf::LogWarning ("CONTROL")
          << "Warning: ParticleGun have encountered larger number of iteration for momentum generation ("
          << iteration << " / " << _iterationLimit << "\n";
      }
    }
    return p;
  }

  double ParticleGunImpl::getHistogramMomentum(void) {
    double urand = _randFlat.fire();
    size_t ibin=0;
    while(urand > histogramBins[ibin]){
      ++ibin;
    }
    return histogramRange[0] + (histogramRange[1]-histogramRange[0])*(ibin+0.5)/float(histogramBins.size());
 }

} // namespace mu2e
