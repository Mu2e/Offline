//
// Cosmic ray muon generator using a TH2 as probability distribution
//
//
// Original author Ralf Ehrlich

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
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/EventGenerator/inc/CosmicFromTH2.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/Mu2eUtilities/inc/rm48.hh"
#include "Offline/GeneralUtilities/inc/safeSqrt.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// ROOT includes
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"

// From CLHEP
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"

using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
using CLHEP::RandFlat;
using CLHEP::GeV;

namespace mu2e
{

  CosmicFromTH2::CosmicFromTH2(CLHEP::HepRandomEngine& engine, art::Run& run, const SimpleConfig& config)
    : _energy(config.getDouble("cosmicFromTH2.energy"))
    , _time(config.getDouble("cosmicFromTH2.time"))
    , _cosmicReferencePointInMu2e (config.getHep3Vector("cosmicFromTH2.cosmicReferencePointInMu2e"))
    , _dx(config.getDouble("cosmicFromTH2.dx"))
    , _dy(config.getDouble("cosmicFromTH2.dy"))
    , _dz(config.getDouble("cosmicFromTH2.dz"))
    , _randFlat{engine}
  {
    GlobalConstantsHandle<ParticleDataList> pdt;
    auto mu_data = pdt->particle(PDGCode::mu_minus);
    double mMu = mu_data.mass();
    _p = safeSqrt(_energy*_energy-mMu*mMu);

    mf::LogInfo log("COSMIC");
    log << "cosmicFromTH2.cosmicReferencePointInMu2e = " << _cosmicReferencePointInMu2e << "\n";
    log << "cosmicFromTH2.histogram = " << _histogram << "\n";
    log << "cosmicFromTH2.energy = " << _energy << "\n";
    log << "cosmicFromTH2.time = " << _time << "\n";

    std::string histogramName = config.getString("cosmicFromTH2.histogram");
    ConfigFileLookupPolicy configFile;
    histogramName = configFile(histogramName);

    TDirectory *directory = gDirectory;
    _file = TFile::Open(histogramName.c_str());
    _histogram = dynamic_cast<TH2*>(gDirectory->FindObjectAny("cosmic"));
    if(_histogram==NULL) throw cet::exception("Configuration")<<"No TH2 histogram found in file.\n";
    directory->cd();

    int n=0;
    if(_dx!=0) n++;
    if(_dy!=0) n++;
    if(_dz!=0) n++;

    if(n!=2) throw cet::exception("Configuration")<<"Production planes not configured correctly. Exactly two out of the dx,dy,dz must be non-zero.\n";
    log << "halflengths: "
        << "cosmicFromTH2.dx = " << _dx <<", "
        << "cosmicFromTH2.dy = " << _dy <<", "
        << "cosmicFromTH2.dz = " << _dz <<"\n";

    // Initialize fake RM48 that is used by DYB code.
    setRm48Distribution(_randFlat);

    // initialize random number generator used by ROOT (for the TH2::GetRandom2 function)
    gRandom->SetSeed(engine.getSeed());
  }  // CosmicDYB()

  CosmicFromTH2::~CosmicFromTH2()
  {
    _file->Close();
  }

  void CosmicFromTH2::generate( GenParticleCollection& genParts )
  {
    if(!_createdProductionPlane)
    {
      GeomHandle<WorldG4>  worldGeom;

      _createdProductionPlane=true;

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

    double theta, phi;
    _histogram->GetRandom2(theta,phi);

    double ct = cos(theta);
    double st = safeSqrt(1. - ct*ct);
    CLHEP::HepLorentzVector mom(_p*st*cos(phi), _p*ct, _p*st*sin(phi), _energy);

    double x = (1.-2.*_randFlat.fire())*_dx;
    double y = (1.-2.*_randFlat.fire())*_dy;  //_dy is 0 for horizontal production planes
    double z = (1.-2.*_randFlat.fire())*_dz;
    CLHEP::Hep3Vector delta(x, y, z);
    CLHEP::Hep3Vector pos = delta + _cosmicReferencePointInMu2e;

    PDGCode::type pid = (_randFlat.fire() > 0.534 ) ? PDGCode::mu_minus : PDGCode::mu_plus;

    genParts.push_back(GenParticle(pid, GenId::cosmicDYB, pos, mom, _time));
  }

}
