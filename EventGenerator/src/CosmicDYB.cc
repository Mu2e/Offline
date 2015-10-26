//
// Cosmic ray muon generator, uses Daya Bay libraries
//
// $Id: CosmicDYB.cc,v 1.32 2014/03/22 21:40:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/03/22 21:40:44 $
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
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
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

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"

// From CLHEP
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::GeV;

namespace mu2e 
{

  CosmicDYB::CosmicDYB( art::Run& run, const SimpleConfig& config )
  : GeneratorBase()

  , _verbose( config.getInt("cosmicDYB.verbose", 0) )
  , _doHistograms( config.getBool("cosmicDYB.doHistograms", true) )
  , _hStartXZ(NULL)
  , _hStartY(NULL)
  , _hStartPlane(NULL)
  , _hStartE(NULL)

  , _mMu(0) //muon mass

    // Mean multiplicity. If negative, use -_mean as a fixed number
  , _mean      ( config.getDouble("cosmicDYB.mean") )
  , _muEMin    ( config.getDouble("cosmicDYB.muEMin") )   //in MeV
  , _muEMax    ( config.getDouble("cosmicDYB.muEMax") )   //in MeV
  , _muCosThMin( config.getDouble("cosmicDYB.muCosThMin") )
  , _muCosThMax( config.getDouble("cosmicDYB.muCosThMax") )

  , _dx(config.getDouble("cosmicDYB.dx"))
  , _dy(0)
  , _dz(config.getDouble("cosmicDYB.dz"))
  , _y0(0)

    // Dimensions of the 2d working space for hrndg2.
  , _ne ( config.getInt("cosmicDYB.nBinsE") )
  , _nth( config.getInt("cosmicDYB.nBinsTheta") )

    // Time range (in ns) over which to generate events.
  ,_tmin( config.getDouble("cosmicDYB.tMin", NAN) )
  ,_tmax( config.getDouble("cosmicDYB.tMax", NAN) )
  ,_dt  ( 0.0 )

    // Random number distributions; getEngine comes from the base class.
  ,_randFlat( getEngine() )
  ,_randPoissonQ( getEngine(), std::abs(_mean) )

    // Working space for hrndg2 (working space will be on the heap).
  , _workingSpace( )
  , _createdProductionPlane(false)

  , _choice(UNDEFINED)
  , _directionChoice(ALL)
  , _cosmicReferencePointInMu2e()
  , _vertical(false)
  , _dontProjectToSurface(config.getBool("cosmicDYB.dontProjectToSurface",false))
  {
    mf::LogInfo log("COSMIC");

    //pick up particle mass
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& mu_data = pdt->particle(PDGCode::mu_minus).ref();
    _mMu = mu_data.mass().value();

    //set _choice to desired cosmic ray generator coordinates
    const std::string refPointChoice = config.getString("cosmicDYB.refPointChoice");
    log << "cosmicDYB.refPointChoice = " << refPointChoice << "\n";

    if (refPointChoice == "Tracker")         _choice = TRACKER;
    else if(refPointChoice == "ExtMonFNAL")  _choice = EXTMONFNAL;
    else if(refPointChoice == "Calorimeter") _choice = CALO;
    else if(refPointChoice == "Customized") 
    {
      _choice = CUSTOMIZED;
      _cosmicReferencePointInMu2e = config.getHep3Vector("cosmicDYB.cosmicReferencePointInMu2e");
      const std::string directionChoice = config.getString("cosmicDYB.directionChoice");
      log << "cosmicDYB.directionChoice = " << directionChoice << "\n";
      if(directionChoice=="All") _directionChoice=ALL;
      else if(directionChoice=="Positive_x") _directionChoice=POSITIVE_X;
      else if(directionChoice=="Negative_x") _directionChoice=NEGATIVE_X;
      else if(directionChoice=="Positive_z") _directionChoice=POSITIVE_Z;
      else if(directionChoice=="Negative_z") _directionChoice=NEGATIVE_Z;
      else throw cet::exception("Configuration")<<"Unknown cosmicDYB.directionChoice\n";
    }
    else 
    {
      throw cet::exception("Configuration") << "Unknown CosmicDYB.refPointChoice\n";
    }

    _vertical = config.getBool("cosmicDYB.vertical", false);
    if(_vertical)
    {
      if(_choice!=CUSTOMIZED) throw cet::exception("Configuration")<<"Vertical production planes require cosmicDYB.refPointChoice: Customized\n";
      if(_directionChoice==ALL) throw cet::exception("Configuration")<<"Vertical production planes require a cosmicDYB.directionChoice other than All e.g. Positive_z or Negative_x\n";
      if(_dx!=0 && _dz!=0) throw cet::exception("Configuration")<<"Vertical production planes must have either cosmicDYB.dx:0 or cosmicDYB.dz:0 \n";
      if(_dx==0 && _dz==0) throw cet::exception("Configuration")<<"Vertical production planes cannot have both cosmicDYB.dx:0 and cosmicDYB.dz:0 \n";
      if((_directionChoice==POSITIVE_X || _directionChoice==NEGATIVE_X) && _dx!=0) throw cet::exception("Configuration")<<"Orientation of the production plane doesn't match cosmicDYB.directionChoice\n";
      if((_directionChoice==POSITIVE_Z || _directionChoice==NEGATIVE_Z) && _dz!=0) throw cet::exception("Configuration")<<"Orientation of the production plane doesn't match cosmicDYB.directionChoice\n";
      _dy=config.getDouble("cosmicDYB.dy");
    }

    if(_choice!=CUSTOMIZED) _y0=config.getDouble("cosmicDYB.y0");

    // Allocate hrndg2 working space on the heap.
    _workingSpace.resize(_ne*_nth);

    log << "cosmicDYB.mean = " << _mean << "\n"
        << "cosmicDYB.muEMin = " << _muEMin <<" MeV, "
        << "cosmicDYB.muEMax = " << _muEMax << " MeV\n"
        << "cosmicDYB.muCosThMin = " << _muCosThMin << ", "
        << "cosmicDYB.muCosThMax = " << _muCosThMax << "\n"
        << "working space dimenions (" << _ne << "," << _nth << ")" << "\n"
        << "cosmicDYB.vertical = " << _vertical <<"\n"
        << "halflengths: "
        << "cosmicDYB.dx = " << _dx <<", "
        << "cosmicDYB.dy = " << _dy <<", "
        << "cosmicDYB.dz = " << _dz <<"\n";
    if(_choice!=CUSTOMIZED) log << "cosmicDYB.y0 = " << _y0 <<"\n";
    else log << "cosmicDYB.cosmicReferencePointInMu2e = " << _cosmicReferencePointInMu2e << "\n";

    // Access conditions data.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");

    // Start time for generation is a little before the start time
    // of the DAQ system.
    double offset = 100.;

    // Start and end times for generation.
    if(isnan(_tmin)) _tmin = (daqPar->t0 > offset)? daqPar->t0-offset : 0.;
    if(isnan(_tmax)) _tmax = accPar->deBuncherPeriod;
    _dt   = _tmax - _tmin;


    // Initialize fake RM48 that is used by DYB code.
    setRm48Distribution(_randFlat);

    // initialize DYB generator
    float par = 1.;

    // convert to GeV
    _muEMin /= GeV;
    _muEMax /= GeV;

    double dim_sum=0;
    double E=0;
    double cosTh=0;
    hrndg2(_workingSpace,_ne,_muEMin,_muEMax,_nth,_muCosThMin,_muCosThMax, dim_sum,E,cosTh,par,_vertical);  //energy is in GeV

    // dim_sum (as returned from the Daya Bay code) is
    // - for horizontal planes: the integral of I(theta,E)*cos(theta)*sin(theta)*dtheta*dE from E0 to E1 and 0 to pi/2 
    // - for vertical planes: the integral of I(theta,E)*sin(theta)^2*dtheta*dE from E0 to E1 and 0 to pi/2 
    // dim_sum has the unit s^-1 * cm^-2
    // I(theta,E) is the intensity from the modified Gaisser formula
    //
    // the rate is 
    // - for horizontal planes: 2*pi*area*integral
    // - for vertical planes: 2*area*integral (only tracks from one direction through the plane are considered)
    //
    // the area is 
    // - for horizontal planes: 4*_dx*_dz
    // - for vertical planes: 4*_dx*_dy or 4*_dz*_dy 
    //
    // the unit is in mm^2 (needs to be multiplied by 0.01 to get it to cm^2 [like in dim_sum])
    // 
    double tRate = dim_sum*M_PI*0.08*_dx*_dz;
    if(!_vertical && _directionChoice!=ALL) tRate = dim_sum*M_PI*0.04*_dx*_dz; // (half the value, since only half of the azimuth angles are considered)
    if(_vertical && _dx!=0) tRate = dim_sum*0.08*_dx*_dy;
    if(_vertical && _dz!=0) tRate = dim_sum*0.08*_dz*_dy;
    log << "Total cosmic rate = " << tRate << " Hz\n";

    if(_directionChoice!=ALL)
    log << "NOTE: The rate above takes into account that only half of the azimuth angles are considered.\n";

    if(_doHistograms)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("CosmicDYB");
      _hStartXZ    = tfdir.make<TH2D>( "StartXZ",    "StartXZ",     500, -6.0e5,  6.0e5, 500, -6.0e5, 6.0e5 );
      _hStartY     = tfdir.make<TH1D>( "StartY",     "StartY",     2000, -5.0e3, 15.0e3 );
      _hStartPlane = tfdir.make<TH1D>( "StartPlane", "StartPlane",    5,  0,      5);
      _hStartE     = tfdir.make<TH1D>( "StartE",     "StartE",      500,  0,      _muEMax*GeV );
      _hStartTheta = tfdir.make<TH1D>( "StartTheta", "StartTheta",  100, -M_PI,   M_PI);
      _hStartPhi   = tfdir.make<TH1D>( "StartPhi",   "StartPhi",    100, -M_PI,   M_PI);
    }
  }  // CosmicDYB()

  CosmicDYB::~CosmicDYB() { }

  void CosmicDYB::generate( GenParticleCollection& genParts )
  {
    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<WorldG4>  worldGeom;

    if(!_createdProductionPlane)
    {
      _createdProductionPlane=true;
    
      GeomHandle<ExtMonFNAL::ExtMon> extMonFNAL;
      GeomHandle<DetectorSystem> detsys;

      switch (_choice)
      {
        case TRACKER:    _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(), _y0, detsys->getOrigin().z());
                         break;
        case EXTMONFNAL: _cosmicReferencePointInMu2e = Hep3Vector(extMonFNAL->detectorCenterInMu2e().x(), _y0, extMonFNAL->detectorCenterInMu2e().z());
                         break;
        case CALO:       _cosmicReferencePointInMu2e = Hep3Vector(detsys->getOrigin().x(), _y0, detsys->getOrigin().z() + 2500.);
                         // distance from tracker to calo center is hardcoded, FIXME!!!!!
                         break;
        case CUSTOMIZED: break;  //already set above
        default:         throw cet::exception("Configuration")<< "Should never occur: unknown CosmicDYB.refPointChoice\n";
                         break;
      }

      const Hep3Vector halfLengths(_dx,_dy,_dz);

      std::cout<<"center of production plane in Mu2e coordinated = "<<_cosmicReferencePointInMu2e<<std::endl;
      std::cout<<"production plane half lengths = "<<halfLengths<<std::endl;
      std::cout<<"Mu2e Origin in the the GEANT world = "<<worldGeom->mu2eOriginInWorld()<<std::endl;
      std::cout<<"GEANT world half lengths = ("
               <<worldGeom->halfLengths()[0]<<", "
               <<worldGeom->halfLengths()[1]<<", "
               <<worldGeom->halfLengths()[2]<<")"<<std::endl;

      if(!worldGeom->inWorld(_cosmicReferencePointInMu2e+halfLengths) || 
         !worldGeom->inWorld(_cosmicReferencePointInMu2e-halfLengths))
      {
        throw cet::exception("GEOM")<<"Cosmic ray production plane is outside of the world volume! Increase the world margins or change production plane\n";
      }
    }

    // Choose the number of muons to generate this event.
    long n = (_mean < 0) ? static_cast<long>(-_mean) : _randPoissonQ.fire();

    for(int i=0; i<n; i++)
    {
      float par = 111.;  // double precision
      double dim_sum,E,cosTh;
      hrndg2(_workingSpace,_ne,_muEMin,_muEMax,_nth,_muCosThMin,_muCosThMax, dim_sum,E,cosTh,par,_vertical);

      // energy is in GeV, convert to MeV
      E *= GeV;

      double p = safeSqrt(E*E-_mMu*_mMu);
      if(E<=_mMu) E = _mMu;

      // Cosine and sin of polar angle wrt y axis.
      double cy = cosTh;
      double sy = safeSqrt(1. - cosTh*cosTh);

      double phi = 2.*M_PI*_randFlat.fire();
      if(_vertical && _dz==0) phi = acos(1-2*_randFlat.fire());             //0...pi
      if(_vertical && _dx==0) phi = acos(1-2*_randFlat.fire()) - M_PI/2.0;  //-pi/2...pi/2

      CLHEP::HepLorentzVector mom(p*sy*cos(phi), -p*cy, p*sy*sin(phi), E);

      switch (_directionChoice)
      {
        case ALL:        break;
        case POSITIVE_X: if(mom.x()<0) mom.setX(-mom.x()); break;
        case NEGATIVE_X: if(mom.x()>0) mom.setX(-mom.x()); break;
        case POSITIVE_Z: if(mom.z()<0) mom.setZ(-mom.z()); break;
        case NEGATIVE_Z: if(mom.z()>0) mom.setZ(-mom.z()); break;
      }

      // Position in reference plane
      double x = (1.-2.*_randFlat.fire())*_dx;
      double y = (1.-2.*_randFlat.fire())*_dy;  //_dy is 0 for horizontal production planes
      double z = (1.-2.*_randFlat.fire())*_dz;
      CLHEP::Hep3Vector delta(x, y, z);
      CLHEP::Hep3Vector pos = delta + _cosmicReferencePointInMu2e;
      if(_verbose>1) std::cout << "position on production plane = " << pos << std::endl;

// project start position (pos) to the surface
      if(!_dontProjectToSurface)
      {
        std::array<double,3> surfaces;

        //surface in y (above the dirt)
        surfaces[1] = env->ymax();

        //surfaces in x and z (world borders in x and z direction on the side in which the track needs to be projected)
        //need the positive x (or z) side if the x (or z) momentum is negative
        CLHEP::Hep3Vector momdir= mom.vect().unit();
        surfaces[0] = (momdir.x()>0?-worldGeom->halfLengths()[0]:worldGeom->halfLengths()[0]) - worldGeom->mu2eOriginInWorld().x();
        surfaces[2] = (momdir.z()>0?-worldGeom->halfLengths()[2]:worldGeom->halfLengths()[2]) - worldGeom->mu2eOriginInWorld().z();

        //find out which surface is hit first
        //for all coordinates i:
        //x_i(plane) = x_i(surface) + t * v_i   where v_i can be substituded by the momentum component
        //and therefore v_i / (x_i(plane) - x_i(surface)) = 1 / t
        //this means, the larger the fraction v_i / (x_i(plane)-x_i(surface)), the larger is 1/t, and the faster the surface i is hit
        //find which coordinate i has the largest fraction
        std::array<double,3> fractions;
        for(int i=0; i<3; i++) fractions[i]=momdir[i]/(pos[i]-surfaces[i]);

        int coord = std::distance(fractions.begin(), std::max_element(fractions.begin(), fractions.end())); //coordinate of largest fraction
        double scale = (surfaces[coord]-pos[coord])/momdir[coord];

        //start position projected to one of the surfaces
        pos += scale*momdir;

        //due to rounding issues, the starting point may have been projected slighlty outside of the world volume
        //need to adjust for this by simply setting this coordinate to the right value
        pos[coord]=surfaces[coord];

        if(_doHistograms)
        {
          switch(coord)
          {
            case 0: _hStartPlane->Fill((momdir.x()>0?"-X":"+X"),1); break;
            case 1: _hStartPlane->Fill((momdir.y()>0?"-Y":"+Y"),1); break;
            case 2: _hStartPlane->Fill((momdir.z()>0?"-Z":"+Z"),1); break;
          }
        }
      }

      if(_doHistograms)
      {
        _hStartXZ->Fill(pos.x(), pos.z());
        _hStartY->Fill(pos.y());
        _hStartE->Fill(E);
        _hStartTheta->Fill(acos(mom.vect().unit().y()));
        _hStartPhi->Fill(atan2(mom.vect().z(),mom.vect().x()));
      }

      if(_verbose>1) std::cout << "starting position = " << pos << std::endl;

      // pick a random starting time, unless a constant time was set
      double time = _tmin + _dt*_randFlat.fire();

      // Pick a random charge.
      // implement a rough charge asymmetry
      double logP = log10(p)-3.;
      double asym = 1.15;
      if(logP>0) 
      {
        if(logP>1) asym = 1.3;
        else asym = 1.15+0.15*logP;
      }

      PDGCode::type pid = (_randFlat.fire() > asym/(1+asym) ) ? PDGCode::mu_minus : PDGCode::mu_plus;

      // Add the muon to the list.
      genParts.push_back(GenParticle(pid, GenId::cosmicDYB, pos, mom, time));

    }
  }

}
