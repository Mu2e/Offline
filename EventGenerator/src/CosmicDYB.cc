//
// Cosmic ray muon generator, uses Daya Bay libraries
//
//

// C++ includes.
#include <cmath>
#include <iostream>

// Framework includes.
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/CosmicDYB.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

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

  CosmicDYB::CosmicDYB(CLHEP::HepRandomEngine& engine, art::Run& run, const SimpleConfig& config )
    : _verbose( config.getInt("cosmicDYB.verbose", 0) )
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
  , _muThMin   ( config.getDouble("cosmicDYB.muThMin", 0.0) )
  , _muThMax   ( config.getDouble("cosmicDYB.muThMax", M_PI/2.0) )
  , _muPhiMin  ( config.getDouble("cosmicDYB.muPhiMin", 0.0) )
  , _muPhiMax  ( config.getDouble("cosmicDYB.muPhiMax", 2.0*M_PI) )

    // Half lengths [mm]
  , _dx(config.getDouble("cosmicDYB.dx"))
  , _dy(config.getDouble("cosmicDYB.dy"))
  , _dz(config.getDouble("cosmicDYB.dz"))

    // Number of lookup bins for DYBGenerator.
  , _ne ( config.getInt("cosmicDYB.nBinsE", 2000) )
  , _nth( config.getInt("cosmicDYB.nBinsTheta", 200) )

    , _randFlat{engine}
    , _randPoissonQ{engine, std::abs(_mean)}

    /*
      TRACKER:(-3904,0,10200)
      EXTMONFNAL:(-642.571,0,-27524.8)
      CALO:(-3904,0,11976.1)
    */
  , _productionCenterInMu2e(config.getHep3Vector("cosmicDYB.productionCenterInMu2e"))
  , _direction(DYBGenerator::Direction::UNDEFINED)
  , _dontProjectToSurface(config.getBool("cosmicDYB.dontProjectToSurface",false))
  , _checkedProductionPlanes(false)
  {
    mf::LogInfo log("COSMIC");

    //pick up particle mass
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const HepPDT::ParticleData& mu_data = pdt->particle(PDGCode::mu_minus).ref();
    _mMu = mu_data.mass().value();


    bool box=false;
    std::string directionString = "Undefined";
    if(_dx!=0 && _dy!=0 && _dz!=0)
    {
       log << "Since dx, dy, and dz are non-zero, the generator uses 5 production planes arranged in a box around the reference point\n";
       log << "The cosmicDYB.direction gets set automatically for all 5 production planes\n";
       log << "The setting for the phi limits get ignored and get set to 0 ... 2*pi\n";
       box=true;
    }
    else
    {
      directionString = config.getString("cosmicDYB.direction");
      if(directionString=="Negative_y")      {_direction=DYBGenerator::Direction::NEGATIVE_Y; _muPhiMin=0.0;       _muPhiMax=2.0*M_PI;}
      else if(directionString=="Positive_x") {_direction=DYBGenerator::Direction::POSITIVE_X; _muPhiMin=-M_PI/2.0; _muPhiMax=M_PI/2.0;}
      else if(directionString=="Positive_z") {_direction=DYBGenerator::Direction::POSITIVE_Z; _muPhiMin=0.0;       _muPhiMax=M_PI;}
      else if(directionString=="Negative_x") {_direction=DYBGenerator::Direction::NEGATIVE_X; _muPhiMin=M_PI/2.0;  _muPhiMax=3*M_PI/2.0;}
      else if(directionString=="Negative_z") {_direction=DYBGenerator::Direction::NEGATIVE_Z; _muPhiMin=M_PI;      _muPhiMax=2.0*M_PI;}
      else throw cet::exception("Configuration")<<"Unknown cosmicDYB.directionChoice\n";

      if(_dx!=0 && _dz!=0)
      {
        if(_direction!=DYBGenerator::Direction::NEGATIVE_Y)
        cet::exception("Configuration")<<"dx!=0 and dz!=0 requires a negative_y direction\n";
      }
      if(_dx!=0 && _dy!=0)
      {
        if(_direction!=DYBGenerator::Direction::NEGATIVE_Z && _direction!=DYBGenerator::Direction::POSITIVE_Z)
        cet::exception("Configuration")<<"dx!=0 and dy!=0 requires a negative_z or positive_z direction\n";
      }
      if(_dy!=0 && _dz!=0)
      {
        if(_direction!=DYBGenerator::Direction::NEGATIVE_X && _direction!=DYBGenerator::Direction::POSITIVE_X)
        cet::exception("Configuration")<<"dy!=0 and dz!=0 requires a negative_x or positive_x direction\n";
      }

      if((_dx==0 && _dy==0) || (_dx==0 && _dz==0) || (_dy==0 && _dz==9))
        cet::exception("Configuration")<<"At least two of the dx, dy, and dz needs to be non-zero\n";

      //azimuth angle can be overriden
      _muPhiMin=config.getDouble("cosmicDYB.muPhiMin", _muPhiMin);
      _muPhiMax=config.getDouble("cosmicDYB.muPhiMax", _muPhiMax);
    }


    log << "cosmicDYB.mean = " << _mean << "\n"
        << "cosmicDYB.muEMin = " << _muEMin <<" MeV, "
        << "cosmicDYB.muEMax = " << _muEMax << " MeV\n"
        << "cosmicDYB.ThMin = " << _muThMin << ", "
        << "cosmicDYB.ThMax = " << _muThMax << "\n"
        << "cosmicDYB.PhiMin = " << _muPhiMin << ", "
        << "cosmicDYB.PhiMax = " << _muPhiMax << "\n"
        << "number of lookup bins (E:" << _ne << ", theta:" << _nth << ")" << "\n"
        << "cosmicDYB.direction = " << directionString << "\n"
        << "halflengths: "
        << "cosmicDYB.dx = " << _dx <<"mm, "
        << "cosmicDYB.dy = " << _dy <<"mm, "
        << "cosmicDYB.dz = " << _dz <<"mm\n"
        << "cosmicDYB.productionCenterInMu2e = " << _productionCenterInMu2e << "\n";

    // Access conditions data.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");

    // convert to GeV
    _muEMin /= GeV;
    _muEMax /= GeV;

    if(box)
    {
      _generators.emplace_back(boost::shared_ptr<DYBGenerator>(new DYBGenerator(DYBGenerator::Direction::NEGATIVE_Y, _muThMin, _muThMax, _muEMin, _muEMax, 0.0,       2.0*M_PI,     _nth, _ne)));
      _generators.emplace_back(boost::shared_ptr<DYBGenerator>(new DYBGenerator(DYBGenerator::Direction::POSITIVE_X, _muThMin, _muThMax, _muEMin, _muEMax, -M_PI/2.0, M_PI/2.0,     _nth, _ne)));
      _generators.emplace_back(boost::shared_ptr<DYBGenerator>(new DYBGenerator(DYBGenerator::Direction::POSITIVE_Z, _muThMin, _muThMax, _muEMin, _muEMax, 0.0,       M_PI,         _nth, _ne)));
      _generators.emplace_back(boost::shared_ptr<DYBGenerator>(new DYBGenerator(DYBGenerator::Direction::NEGATIVE_X, _muThMin, _muThMax, _muEMin, _muEMax, M_PI/2.0,  3.0/2.0*M_PI, _nth, _ne)));
      _generators.emplace_back(boost::shared_ptr<DYBGenerator>(new DYBGenerator(DYBGenerator::Direction::NEGATIVE_Z, _muThMin, _muThMax, _muEMin, _muEMax, M_PI,      2.0*M_PI,     _nth, _ne)));
    }
    else
    {
      _generators.emplace_back(boost::shared_ptr<DYBGenerator>(new DYBGenerator(_direction, _muThMin, _muThMax, _muEMin, _muEMax, _muPhiMin, _muPhiMax, _nth, _ne)));
    }

    double totalRate=0;
    for(unsigned int i=0; i<_generators.size(); i++)
    {
      double rate = _generators[i]->GetRate();  //in cm^-2 * s^-1
      switch(_generators[i]->GetDirection())
      {
        case DYBGenerator::Direction::NEGATIVE_Y: rate*=0.04*_dx*_dz; break;
        case DYBGenerator::Direction::POSITIVE_X:
        case DYBGenerator::Direction::NEGATIVE_X: rate*=0.04*_dy*_dz; break;
        case DYBGenerator::Direction::POSITIVE_Z:
        case DYBGenerator::Direction::NEGATIVE_Z: rate*=0.04*_dx*_dy; break;
                                         default: throw cet::exception("Configuration")<<"Invalid direction in DYB generator\n";
      };
      _boxFraction.push_back(rate);
      totalRate+=rate;
    }
    for(unsigned int i=0; i<_generators.size(); i++)
    {
      _boxFraction[i]/=totalRate;
      if(i>0) _boxFraction[i]+=_boxFraction[i-1];
    }
    _boxFraction.back()=1; //to counteract rounding errors;

    log << "Total cosmic rate = " << totalRate << " Hz\n";

    if(_doHistograms)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("CosmicDYB");
      _hStartXZ    = tfdir.make<TH2D>( "StartXZ",    "StartXZ",    1000, -6.0e5,  6.0e5, 1000, -6.0e5, 6.0e5 );
      _hStartY     = tfdir.make<TH1D>( "StartY",     "StartY",     1000, -5.0e3, 20.0e3 );
      _hStartPlane = tfdir.make<TH1D>( "StartPlane", "StartPlane",    5,  0,      5);
      _hStartE     = tfdir.make<TH1D>( "StartE",     "StartE",     1000,  0,      1.0e6);
      _hStartTheta = tfdir.make<TH1D>( "StartTheta", "StartTheta", 1000, -M_PI,   M_PI);
      _hStartPhi   = tfdir.make<TH1D>( "StartPhi",   "StartPhi",   1000, -M_PI,   M_PI);
    }

  }  // CosmicDYB()

  CosmicDYB::~CosmicDYB() { }

  void CosmicDYB::generate( GenParticleCollection& genParts )
  {
    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<WorldG4>  worldGeom;

    if(!_checkedProductionPlanes)
    {
      _checkedProductionPlanes=true;

      const Hep3Vector halfLengths(_dx,_dy,_dz);

      std::cout<<"production center in Mu2e coordinates = "<<_productionCenterInMu2e<<std::endl;
      std::cout<<"production plane half lengths = "<<halfLengths<<std::endl;
      std::cout<<"Mu2e Origin in the the GEANT world = "<<worldGeom->mu2eOriginInWorld()<<std::endl;
      std::cout<<"GEANT world half lengths = ("
               <<worldGeom->halfLengths()[0]<<", "
               <<worldGeom->halfLengths()[1]<<", "
               <<worldGeom->halfLengths()[2]<<")"<<std::endl;

      if(!worldGeom->inWorld(_productionCenterInMu2e+halfLengths) ||
         !worldGeom->inWorld(_productionCenterInMu2e-halfLengths))
      {
        throw cet::exception("GEOM")<<"Cosmic ray production plane is outside of the world volume! Increase the world margins or change production plane\n";
      }
    }

    // Choose the number of muons to generate this event.
    long n = (_mean < 0) ? static_cast<long>(-_mean) : _randPoissonQ.fire();

    for(int i=0; i<n; i++)
    {
      double theta, phi, E;
      unsigned int side=0;
      if(_generators.size()==1) _generators[0]->GenerateMuon(theta, E, phi, _randFlat);
      else
      {
        double uSide = _randFlat.fire();
        for(; side<_boxFraction.size(); side++)
        {
          if(uSide<_boxFraction[side]) {_generators[side]->GenerateMuon(theta, E, phi, _randFlat); break;}
        }
      }

      // energy is in GeV, convert to MeV
      E *= GeV;

      double p = safeSqrt(E*E-_mMu*_mMu);   //TODO: Is E really the total energy?
      if(E<=_mMu) E = _mMu;

      // Cosine and sin of polar angle wrt y axis.
      double cosTh = cos(theta);
      double sinTh = safeSqrt(1. - cosTh*cosTh);

      CLHEP::HepLorentzVector mom(p*sinTh*cos(phi), -p*cosTh, p*sinTh*sin(phi), E);

      // Position in reference plane
      double x = (1.-2.*_randFlat.fire())*_dx;
      double y = (1.-2.*_randFlat.fire())*_dy;
      double z = (1.-2.*_randFlat.fire())*_dz;
      if(_generators.size()>1)
      {
        switch(_generators[side]->GetDirection())
        {
          case DYBGenerator::Direction::NEGATIVE_Y: y= _dy; break;
          case DYBGenerator::Direction::POSITIVE_X: x=-_dx; break;
          case DYBGenerator::Direction::NEGATIVE_X: x= _dx; break;
          case DYBGenerator::Direction::POSITIVE_Z: z=-_dz; break;
          case DYBGenerator::Direction::NEGATIVE_Z: z= _dz; break;
                                           default: throw cet::exception("Configuration")<<"Invalid direction in DYB generator\n";
        };
      }

      CLHEP::Hep3Vector delta(x, y, z);
      CLHEP::Hep3Vector pos = delta + _productionCenterInMu2e;
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

      // Add the muon to the list. time = 0.0
      genParts.push_back(GenParticle(pid, GenId::cosmicDYB, pos, mom, 0.0));

    }
  }

}
