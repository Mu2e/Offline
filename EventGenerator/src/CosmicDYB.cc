//
// Cosmic ray muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.cc,v 1.31 2014/03/19 18:47:29 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/03/19 18:47:29 $
//
// Original author Yury Kolomensky
//
// Notes:
// 1) The hrndg2 generator takes a very long time during initialization if the
//    22 seconds on ilcsim2 if muEMin = 10,000 MeV.
//    31 seconds on ilcsim2 if muEMin =  5,000 MeV.
//    41 seconds on ilcsim2 if muEMin =  3,001 MeV.
//    52 seconds on ilcsim2 if muEmin =  1,000 MeV (same for 1,001 MeV and for 999 MeV).
//    68 seconds on ilcsim2 if muEmin =    601 MeV
//
//    45 seconds on ilcsim  if muEMin = 10,001 MeV.
//    62 seconds on ilcsim  if muEMin =  5,001 MeV.
//

// C++ includes.
#include <cmath>
#include <iostream>

// Framework includes.
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "EventGenerator/inc/CosmicDYB.hh"
#include "EventGenerator/inc/hrndg2.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/rm48.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"


// From CLHEP
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

// From Root.
#include "TH1D.h"
#include "TH2D.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::GeV;

namespace mu2e {

  // Mass of the muon, in MeV.
  // Once we have the HepPDT package installed, get the mass from there.
  static const double mMu = 105.6584;

  CosmicDYB::CosmicDYB( art::Run& run, const SimpleConfig& config )
  : GeneratorBase()
    // Histograms
  , _cosmicMultiplicityH( 0 )
  , _cosmicMomH         ( 0 )
  , _cosmicChargeH      ( 0 )
  , _cosmicCosThetaH    ( 0 )
  , _cosmicCosThetaVsEH ( 0 )

    // configurable parameters

    // Mean multiplicity. If negative, use -_mean as a fixed number
  , _mean      ( config.getDouble("cosmicDYB.mean",0.) )
  , _muEMin    ( config.getDouble("cosmicDYB.muEMin", 3000.) )
  , _muEMax    ( config.getDouble("cosmicDYB.muEMax", 1.e5) )
  , _muCosThMin( config.getDouble("cosmicDYB.muCosThMin",0.00366518) )
  , _muCosThMax( config.getDouble("cosmicDYB.muCosThMax",1.0) )

    // Dimensions of the 2d working space for hrndg2.
  , _ne ( config.getInt("cosmicDYB.nBinsE", _default_ne) )
  , _nth( config.getInt("cosmicDYB.nBinsTheta",_default_nth) )

    // Control of histograms.
  , _doHistograms( config.getBool("cosmicDYB.doHistograms", true) )
  , _verbose( config.getInt("cosmicDYB.verbose", 0) )

    // end of configurable parameters

    // Time range (in ns) over which to generate events.
  ,_tmin( 0.0 )
  ,_tmax( 0.0 )
  ,_dt  ( 0.0 )

    // Random number distributions; getEngine comes from the base class.
  ,_randFlat( getEngine() )
  ,_randPoissonQ( getEngine(), std::abs(_mean) )

    // Working space for hrndg2 (working space will be on the heap).
  , _workingSpace( )
  , _checkedProductionPlane(false)


  , _choice(UNDEFINED)
  , _directionChoice(ALL)

  {
    mf::LogInfo log("COSMIC");

    _dx=config.getDouble("cosmicDYB.dx",5000);
    _dz=config.getDouble("cosmicDYB.dz",5000);
    _y0=config.getDouble("cosmicDYB.y0",0);

    // Sanity check.
    if ( std::abs(_mean) > 99999. ) {
      throw cet::exception("RANGE")
        << "CosmicDYB has been asked to produce a crazily large number of electrons."
        << _mean
        << "\n";
    }

    //set _choice to desired cosmic ray generator coordinates
    const std::string refPointChoice = config.getString("cosmicDYB.refPointChoice", "Customized");
    log << "cosmicDYB.refPointChoice = " << refPointChoice << "\n";
    if (refPointChoice == "Tracker") {
      _choice = TRACKER;
    }

    else if(refPointChoice == "ExtMonFNAL") {
      _choice = EXTMONFNAL;
    }

    else if(refPointChoice == "Calorimeter") {
      _choice = CALO;
    }

    else if(refPointChoice == "Customized") {
      _choice = CUSTOMIZED;
      _cosmicReferencePointInMu2e = config.getHep3Vector("cosmicDYB.cosmicReferencePointInMu2e");
      const std::string directionChoice = config.getString("cosmicDYB.directionChoice");
      log << "cosmicDYB.directionChoice = " << directionChoice << "\n";
      if(directionChoice=="ALL") _directionChoice=ALL;
      else if(directionChoice=="Positive_z") _directionChoice=POSITIVE_Z;
      else if(directionChoice=="Negative_x") _directionChoice=NEGATIVE_X;
      else throw cet::exception("Configuration")<<"Unknown cosmicDYB.directionChoice\n";
    }

    else {
      throw cet::exception("Configuration")
        << "Unknown CosmicDYB.refPointChoice"
        << "\n";
    }


    // Allocate hrndg2 working space on the heap.
    _workingSpace.resize(_ne*_nth);

    log << "cosmicDYB.mean = " << _mean << "\n"
        << "cosmicDYB.muEMin = " << _muEMin
        << ", cosmicDYB.muEMax = " << _muEMax << "\n"
        << "cosmicDYB.muCosThMin = " << _muCosThMin
        << ", cosmicDYB.muCosThMax = " << _muCosThMax << "\n"
        << "cosmicDYB.dx = " << _dx
        << ", cosmicDYB.dz = " << _dz
        << ", cosmicDYB.y0 = " << _y0
        << ", working space dimenions ("
        << _ne << "," << _nth << ")" << "\n";

    // Access conditions data.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");

    // Start time for generation is a little before the start time
    // of the DAQ system.
    double offset = 100.;

    // Start and end times for generation.
    _tmin = (daqPar->t0 > offset)? daqPar->t0-offset : 0.;
    _tmax = accPar->deBuncherPeriod;
    _dt   = _tmax - _tmin;


    // Book histograms in a separate subdirectory.
    if ( _doHistograms ){

      art::ServiceHandle<art::TFileService> tfs;

      art::TFileDirectory tfdir = tfs->mkdir( "CosmicDYB" );
      _cosmicMultiplicityH = tfdir.make<TH1D>( "MultiplicityH", "Cosmic Multiplicity", 20, -0.5, 19.5);

      // log of muon energy (GeV)
      _cosmicMomH = tfdir.make<TH1D>( "MomH", "log (Momentum, GeV)", 60, -3, 3);

      // charge
      _cosmicChargeH = tfdir.make<TH1D>( "ChargeH", "Muon Charge", 2, -2, 2.);

      // cos(theta)
      _cosmicCosThetaH = tfdir.make<TH1D>( "CosThetaH", "Cos(Theta)", 60, -1, 1.);

      // cos(theta) vs log(energy)
      _cosmicCosThetaVsEH = tfdir.make<TH2D>( "CosThetavsEH", "Cos(Theta) vs log (Momentum, GeV)",
                                              60, -3., 3., 60, -1., 1.);
    }

    // Initialize fake RM48 that is used by DYB code.
    setRm48Distribution(_randFlat);

    // initialize DYB generator
    float par = 1.;

    // convert to GeV
    _muEMin /= GeV;
    _muEMax /= GeV;

    double dim_sum,E,cosTh;
    hrndg2(_workingSpace,_ne,_muEMin,_muEMax,_nth,_muCosThMin,_muCosThMax,
           dim_sum,E,cosTh,par);

    // rate is per cm^2. The constants are 2*CLHEP::pi times the area
    double tRate = dim_sum*M_PI*0.08*_dx*_dz;
    log << "Total cosmic rate = " << tRate << " Hz\n";
    if(_directionChoice==POSITIVE_Z || _directionChoice==NEGATIVE_X)
    log << "NOTE: Since only half of the azimuth angles are considered, the rate is only half of the above value.\n";

  }  // CosmicDYB()

  CosmicDYB::~CosmicDYB() { }

  void CosmicDYB::generate( GenParticleCollection& genParts )
  {
    GeomHandle<ExtMonFNAL::ExtMon>  extMonFNAL;
    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<DetectorSystem> detsys;

    // here: take different values depending on config
    // e.g. extMonFNAL->detectorCenterInMu2e()
    switch (_choice){
    case TRACKER:
      _cosmicReferencePointInMu2e =  Hep3Vector(detsys->getOrigin().x(),
                                         _y0,
                                         detsys->getOrigin().z());
      break;
    case EXTMONFNAL:
      _cosmicReferencePointInMu2e = Hep3Vector(extMonFNAL->detectorCenterInMu2e().x(),
                                              _y0,
                                              extMonFNAL->detectorCenterInMu2e().z());
      break;
    case CALO:
    // distance from tracker to calo center is hardcoded, FIXME!!!!!
      _cosmicReferencePointInMu2e =  Hep3Vector(detsys->getOrigin().x(),
                                         _y0,
                                         detsys->getOrigin().z() + 2500.);
      break;
    case CUSTOMIZED:
      break;
    default:
      throw cet::exception("Configuration")
        << "Should never occur: unknown CosmicDYB.refPointChoice"
        << "\n";
      break;
    }

    if(!_checkedProductionPlane)
      {
        _checkedProductionPlane=true;

        std::cout<<"cosmicDYB: cosmicReferencePointInMu2e = "<<_cosmicReferencePointInMu2e<<std::endl;

        art::ServiceHandle<GeometryService> geom;

        if(geom->hasElement<WorldG4>()) {

          GeomHandle<WorldG4>  worldGeom;
          Hep3Vector const& mu2eOrigin = worldGeom->mu2eOriginInWorld();
          std::vector<double> const& halfLengths = worldGeom->halfLengths();

          double marginXMin = halfLengths[0] + (_cosmicReferencePointInMu2e.x()-_dx + mu2eOrigin.x());
          double marginZMin = halfLengths[2] + (_cosmicReferencePointInMu2e.z()-_dz + mu2eOrigin.z());
          double marginXMax = halfLengths[0] - (_cosmicReferencePointInMu2e.x()+_dx + mu2eOrigin.x());
          double marginZMax = halfLengths[2] - (_cosmicReferencePointInMu2e.z()+_dz + mu2eOrigin.z());
          double marginYMin = halfLengths[1] + (_cosmicReferencePointInMu2e.y()+env->ymax() + mu2eOrigin.y());
          double marginYMax = halfLengths[1] - (_cosmicReferencePointInMu2e.y()+env->ymax() + mu2eOrigin.y());

          if(_verbose>0)std::cout<<std::endl<<"distances from the edges of the cosmic ray production plane to the borders of the world volume:"<<std::endl
                   <<"in negative x direction: "<<marginXMin<<std::endl
                   <<"in positive x direction: "<<marginXMax<<std::endl
                   <<"in negative y direction: "<<marginYMin<<std::endl
                   <<"in positive y direction: "<<marginYMax<<std::endl
                   <<"in negative z direction: "<<marginZMin<<std::endl
                   <<"in positive z direction: "<<marginZMax<<std::endl<<std::endl;

          checkCosmicRayProductionPlane(marginXMin, "world.margin.xmin");
          checkCosmicRayProductionPlane(marginXMax, "world.margin.xmax");
          checkCosmicRayProductionPlane(marginYMin, "world.margin.top");
          checkCosmicRayProductionPlane(marginYMax, "world.margin.bottom");
          checkCosmicRayProductionPlane(marginZMin, "world.margin.zmin");
          checkCosmicRayProductionPlane(marginZMax, "world.margin.zmax");
        }
      }


    // Choose the number of electrons to generate this event.
    long n = (_mean < 0) ? static_cast<long>(-_mean) : _randPoissonQ.fire();
    if ( _doHistograms) {
      _cosmicMultiplicityH->Fill(n);
    }

    for ( int i = 0; i != n; ++i ){

      float par = 111.;  // double precision
      double dim_sum,E,cosTh;
      hrndg2( _workingSpace,_ne,_muEMin,_muEMax,_nth,_muCosThMin,_muCosThMax,
              dim_sum,E,cosTh,par);

      // energy is in GeV, convert to MeV
      E *= GeV;

      double p = safeSqrt(E*E-mMu*mMu);
      if ( E <= mMu ) {
        E = mMu;
      }

      // Log10 of momentum, in GeV
      double log10p = log10(p)-3.;

      _cosmicMomH->Fill(log10p);
      _cosmicCosThetaH->Fill(cosTh);
      _cosmicCosThetaVsEH->Fill(log10p,cosTh);

      // Cosine and sin of polar angle wrt y axis.
      double cy = cosTh;
      double sy = safeSqrt(1. - cosTh*cosTh);

      double phi = 2.*M_PI*_randFlat.fire();

      CLHEP::HepLorentzVector mom(p*sy*cos(phi), -p*cy, p*sy*sin(phi), E);

      switch (_directionChoice)
      {
        case POSITIVE_Z: 
             if(mom.z()<0) mom.setZ(-mom.z());
             break;
        case NEGATIVE_X: 
             if(mom.x()>0) mom.setX(-mom.x());
             break;
      }

      // Position in reference plane
      double x = (1.-2.*_randFlat.fire())*_dx;
      double y = 0.0;
      double z = (1.-2.*_randFlat.fire())*_dz;
      CLHEP::Hep3Vector pos( x, y, z );
      CLHEP::Hep3Vector refposInMu2eCoordinates=pos + _cosmicReferencePointInMu2e;

// move this to a position just above the ground
      double ymax = env->ymax();
      CLHEP::Hep3Vector momdir= mom.vect().unit();
      double scale = (ymax-refposInMu2eCoordinates.y())/momdir.y();
      CLHEP::Hep3Vector shift = scale*momdir;
      CLHEP::Hep3Vector posInMu2eCoordinates = refposInMu2eCoordinates + shift;
      if(_verbose > 1){
	cout << "position in reference plane = " << refposInMu2eCoordinates << endl;
	cout << "position at origin = " << posInMu2eCoordinates << endl;
      }

      double time = _tmin + _dt*_randFlat.fire();

      // Pick a random charge.
      // implement a rough charge asymmetry
      double logP = log10(p)-3.;
      double asym = 1.15;
      if ( logP > 0 ) {
        if ( logP > 1 ) {
          asym = 1.3;
        } else {
          asym = 1.15+0.15*logP;
        }
      }


      PDGCode::type pid = (_randFlat.fire() > asym/(1+asym) ) ? PDGCode::mu_minus : PDGCode::mu_plus;

      _cosmicChargeH->Fill(-pid/abs(pid));

      // Add the cosmic to  the list.
      if(filterGeneratedMuons(posInMu2eCoordinates))
      genParts.push_back( GenParticle( pid, GenId::cosmicDYB, posInMu2eCoordinates, mom, time));

    }
  }

  void CosmicDYB::checkCosmicRayProductionPlane(double margin, std::string name)
  {
    if(margin<0) throw cet::exception("GEOM")<<"Cosmic ray production plane is outside of the world volume! \n"
                                             <<"Increase "<<name<<" by at least "<<-margin<<"\n";
  }

  bool CosmicDYB::filterGeneratedMuons(CLHEP::Hep3Vector const &posInMu2eCoordinates)
  {

    art::ServiceHandle<GeometryService> geom;

    if(geom->hasElement<WorldG4>()) {

      GeomHandle<WorldG4>  worldGeom;
      Hep3Vector const& mu2eOrigin = worldGeom->mu2eOriginInWorld();
      std::vector<double> const& halfLengths = worldGeom->halfLengths();
      Hep3Vector diff = posInMu2eCoordinates-mu2eOrigin;
      static const double margin(10.0); // 10mm margin
      if(halfLengths[0] - fabs(diff.x()) < margin ||
	  halfLengths[1] - fabs(diff.y()) < margin ||
	  halfLengths[2] - fabs(diff.z()) < margin) return false;
    }

    return true;
  }
}
