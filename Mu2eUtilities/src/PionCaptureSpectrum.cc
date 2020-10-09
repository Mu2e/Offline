// Simple approximations available for DIO spectrum.
//
//

#include <cmath>
#include <stdlib.h>
// C++ includes
#include <iostream>
#include <memory>

// CLHEP includes
#include "CLHEP/Random/RandFlat.h"

#include "CLHEP/Vector/LorentzVector.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "DataProducts/inc/PDGCode.hh"
// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "HepPDT/Measurement.hh"

#include "HepPDT/ParticleData.hh"
#include "cetlib/pow.h"

#include "Mu2eUtilities/inc/PionCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

using namespace std;

using cet::pow;
using cet::square;

namespace mu2e {

  PionCaptureSpectrum::PionCaptureSpectrum(CLHEP::RandFlat* rnFlat, RandomUnitSphere* rnUnitSphere):
    _spectrum    (Bistirlich      ),
    _spectrum2D  (KrollWadaJoseph ),
    _kMaxUserSet (false           ),
    _kMaxUser    (0.              ),
    _kMaxMax     (0.              ),
    _rnFlat      (rnFlat          ),
    _rnUnitSphere(rnUnitSphere    )
  {
    GlobalConstantsHandle<ParticleDataTable> pdt;

    _me    = pdt->particle(PDGCode::e_minus ).ref().mass().value();
    _mpi   = pdt->particle(PDGCode::pi_minus).ref().mass().value();
    _MN    = GlobalConstantsHandle<PhysicsParams>()->getAtomicMass("Al");
  }

  PionCaptureSpectrum::PionCaptureSpectrum(bool kMaxUserSet, double kMaxUser, double kMaxMax,
					   CLHEP::RandFlat* rnFlat, RandomUnitSphere* rnUnitSphere):
    _spectrum    (Bistirlich      ),
    _spectrum2D  (KrollWadaJoseph ),
    _kMaxUserSet (kMaxUserSet     ),
    _kMaxUser    (kMaxUser        ),
    _kMaxMax     (kMaxMax         ),
    _rnFlat      (rnFlat          ),
    _rnUnitSphere(rnUnitSphere    )
  {
    GlobalConstantsHandle<ParticleDataTable> pdt;

    _me    = pdt->particle(PDGCode::e_minus ).ref().mass().value();
    _mpi   = pdt->particle(PDGCode::pi_minus).ref().mass().value();
    _MN    = GlobalConstantsHandle<PhysicsParams>()->getAtomicMass("Al");
  }

  double PionCaptureSpectrum::getWeight(double E) const {

    double weight(0.);
    if      ( _spectrum == Flat )       weight = getFlat( E );
    else if ( _spectrum == Bistirlich ) weight = getBistirlichSpectrum( E );

    return weight;

  }

  double PionCaptureSpectrum::get2DWeight(double x, double y, double E) const {

    double weight(0.);
    if      ( _spectrum2D == Flat2D          ) weight = getFlat( E, x, y );
    else if ( _spectrum2D == KrollWadaJoseph ) weight = getKrollWadaJosephSpectrum( E, x, y );

    return weight;

  }

  //======================================================
  // Flat spectrum to be used with weighting
  //======================================================
  double PionCaptureSpectrum::getFlat(double E, double x, double y) const {
    return 1.;
  }

//=======================================================
// Analytic fit to the photon energy spectrum for Mg
//  - J.A.Bistirlich, et al, Phys. Rev. C 5, 1867 (1972)
//  - commented out: fit results given in I. Sarra, mu2e-docdb-665-v2
//    Ivano fitted the distribution in fig 7b in the paper, whereas the one
//    to fit is the one in fig 7a. The paper is not very clear on that,
//    see mu2e-19366 for details
//    new fit is the same parameterization of the spectrum in fig 7a,
//    normalization the fit function: integral(f) = 1
//=======================================================
  double PionCaptureSpectrum::getBistirlichSpectrum(double e) const {
    // static const double emax  = 138.2;
    // static const double alpha =   2.691;
    // static const double gamma =   1.051;
    // static const double tau   =   8.043;
    // static const double c0    =   2.741;
    // static const double c1    =  -0.005;
    // return std::pow(emax-e,alpha) * std::exp(-(emax-gamma*e)/tau) * (c0 + c1*e) / 10000;

    static const double emax  =   1.34530e2;
    static const double alpha =   1.29931;
    static const double gamma =   9.28705e-1;
    static const double tau   =   9.39676;
    static const double c0    =   4.77611e-2;
    static const double c1    =  -3.26349e-4;

    if ( e > emax ) return 0.;

    return std::pow(emax-e,alpha) * std::exp(-(emax-gamma*e)/tau) * (c0+c1*e);
  }

  //=======================================================
  // Analytic expression for electron/positron energies via
  //  - Kroll and Wada, Phys. Rev. 98, 1355 (1955)
  //  - Joseph, Il Nuovo Cimento 16, 997 (1960)
  //=======================================================

  double PionCaptureSpectrum::get2DMax(double E) const {
    static const GlobalConstantsHandle<ParticleDataTable> pdt;
    static const HepPDT::ParticleData& e_data  = pdt->particle(PDGCode::e_minus).ref();

    return getKrollWadaJosephSpectrum( E, 2*e_data.mass().value(), 0. );
  }


//------------------------------------------------------------------------------
// E - energy (k0) of the emitted virtual photon, x - mass of the e+e- pair
//-----------------------------------------------------------------------------
  double PionCaptureSpectrum::getKrollWadaJosephSpectrum(double E, double x, double y) const {


    // mass of the recoiling system, neglecting the muon binding energy and nuclear recoil

    double M = _MN + _mpi - E;
					// Set pdf to zero if x is out of bounds
    if (  x > E || x < 2*_me ) return 0.;

					// Set pdf to zero if x or y are out of bounds
    double eta = sqrt( 1 - pow(2*_me/x,2) );
    if ( abs( y ) > eta ) return 0.;

    // Set parameters
    static const double rV = 0.57;
    static const double muS = 0.064;

    double xe = x/E;
    double rT = 1+rV*rV/3 - 2*muS*xe*xe;

    double xl = (1./(1-0.466*xe*xe)+1.88*(rV*rV/6 + muS*(1-xe*xe)))/(1+muS);
    double rL = 0.142*xl*xl;

    double kRatio2 = ( pow(2*E*M + E*E,2) - 2*x*x*(2*M*M + 2*E*M + E*E ) + pow( x, 4 ) )/pow(2*E*M + E*E,2) ;
    double kRatio  = sqrt( kRatio2 );

    double prefactor = ( pow( E+M,2 ) + M*M - x*x )/( pow( E+M,2 ) + M*M );

    double trans     = rT*( ( 1+y*y )/x + 4*_me*_me/pow(x,3) );
    double longit    = rL*(1-y*y)*8*pow(E+M,2)*x/pow( 2*E*M + E*E + x*x, 2 );

    double prob      = kRatio*prefactor*( trans + longit );

    return prob;

  }

//-----------------------------------------------------------------------------
// return X and Y with the probability defined by energy
//-----------------------------------------------------------------------------
  void PionCaptureSpectrum::fire(double energy, double& x, double& y) const {

    double pdfMax = get2DMax(energy);
    x = 0;
    y = 0;
    double prob(0.), threshold(0.);
    do {
      x         = 2*_me + ( energy-2*_me)*_rnFlat->fire();
      y         = -1.   + 2.*_rnFlat->fire();
      threshold = get2DWeight(x,y,energy);
      prob      = pdfMax*_rnFlat->fire();
    } while (prob > threshold);
  }

//-----------------------------------------------------------------------------
  void PionCaptureSpectrum::getElecPosiVectors(double ePhoton, CLHEP::HepLorentzVector& mome, CLHEP::HepLorentzVector& momp) const {

    double x, y;
					// generate invariant mass and energy splitting
    fire(ePhoton,x,y);

    // Get electron/positron energies from x, y values (see Mu2eUtilities/src/PionCaptureSpectrum.cc for details)
    double eElectron = 0.5*( ePhoton - y*std::sqrt( cet::diff_of_squares( ePhoton, x ) ) );
    double ePositron = 0.5*( ePhoton + y*std::sqrt( cet::diff_of_squares( ePhoton, x ) ) );

    // Get electron/positron momentum magnitudes
    double pElectron = std::sqrt( cet::diff_of_squares(eElectron, _me) );
    double pPositron = std::sqrt( cet::diff_of_squares(ePositron, _me) );

    // Produce electron momentum
    CLHEP::Hep3Vector p3electron = _rnUnitSphere->fire( pElectron );

    // Get positron momentum
    CLHEP::Hep3Vector p3positron( p3electron );
    p3positron.setMag( pPositron );

    // - theta (opening angle wrt electron) is constrained by virtual mass formula
    // - phi is allowed to vary between 0 and 2pi

    double cosTheta = 1/(2*pElectron*pPositron)*( cet::square(ePhoton) - cet::sum_of_squares( x, pElectron, pPositron) );

    double phi = 2*M_PI*_rnFlat->fire();

    // - find a vector that is not collinear with the electron direction
    CLHEP::Hep3Vector n1 = (std::abs(p3electron.x()) < std::abs(p3electron.y())) ?
      ((std::abs(p3electron.x()) < std::abs(p3electron.z())) ? CLHEP::Hep3Vector(1,0,0) : CLHEP::Hep3Vector(0,0,1)) :
      ((std::abs(p3electron.x()) < std::abs(p3electron.y())) ? CLHEP::Hep3Vector(1,0,0) : CLHEP::Hep3Vector(0,1,0));

    // - construct a vector perpendicular to the electron momentum
    CLHEP::Hep3Vector perp = p3electron.cross(n1);

    p3positron.rotate(perp      , std::acos(cosTheta) );
    p3positron.rotate(p3electron, phi                 );

    mome.set(p3electron,eElectron);
    momp.set(p3positron,ePositron);
  }

}
