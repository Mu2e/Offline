// C++ includes.
#include <cmath>
#include <iostream>
#include <cstdlib>

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
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/hrndg2.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/rm48.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "CRYGenerator.h"
#include "CRYSetup.h"
#include "CRYEventGenerator/inc/CosmicCRY.hh"
// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

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
   // The following part is needed for the RNG
  template<class T> class RNGWrapper {
    public:
      static void set(T* object, double (T::*func)(void));
      static double rng(void);
    private:
      static T* m_obj;
      static double (T::*m_func)(void);
  };// end of RNGWrapper class

  template<class T> T* RNGWrapper<T>::m_obj;
  template<class T> double (T::*RNGWrapper<T>::m_func)(void);
  template<class T> void RNGWrapper<T>::set(T* object, double (T::*func)(void)) {
    m_obj = object; m_func = func;
  }
  template<class T> double RNGWrapper<T>::rng(void) { return (m_obj->*m_func)(); }
            

  CosmicCRY::CosmicCRY( art::Run& run, const SimpleConfig& config )
  : _verbose(config.getInt("cosmicCRY.verbose", 0) )
  , _doHistograms(config.getBool("cosmicCRY.doHistograms", true) )
  , _saveTree(config.getBool("cosmicCRY.saveTree", false) )
  , _hStartXZ(NULL)
  , _hStartY(NULL)
  // , _hStartPlane(NULL)
  , _hStartE(NULL)
  , _hStartTheta(NULL)
  , _hStartPhi(NULL)

  , _muEMin    ( config.getDouble("cosmicCRY.muEMin") )   //in MeV
  , _muEMax    ( config.getDouble("cosmicCRY.muEMax") )   //in MeV
  , _muCosThMin( config.getDouble("cosmicCRY.muCosThMin") )
  , _muCosThMax( config.getDouble("cosmicCRY.muCosThMax") )
  , _muPhiMin(0)
  , _muPhiMax(2.0*M_PI)

  , _returnMuons( config.getBool("cosmicCRY.returnMuons", true))
  , _returnNeutrons( config.getBool("cosmicCRY.returnNeutrons", true))
  , _returnProtons( config.getBool("cosmicCRY.returnProtons", true))
  , _returnGammas( config.getBool("cosmicCRY.returnGammas", true))
  , _returnElectrons( config.getBool("cosmicCRY.returnElectrons", true))
  , _returnPions( config.getBool("cosmicCRY.returnPions", true))
  , _returnKaons( config.getBool("cosmicCRY.returnKaons", true))

  , _month(config.getInt("cosmicCRY.month", 6))
  , _day(config.getInt("cosmicCRY.day", 21))
  , _year(config.getInt("cosmicCRY.year", 2020))
  , _latitude( config.getDouble("cosmicCRY.latitude", 41.8))
  , _altitude( config.getInt("cosmicCRY.altitude", 0))
  , _subboxLength( config.getDouble("cosmicCRY.subboxLength", 100.))

  , _setupString("")

  {
    mf::LogInfo log("CRYCOSMIC");

    if (_returnMuons) 
      _setupString.append("returnMuons 1 "); // must have the trailing white space
    else
      _setupString.append("returnMuons 0 ");

    if (_returnNeutrons) 
      _setupString.append("returnNeutrons 1 ");
    else
      _setupString.append("returnNeutrons 0 ");

    if (_returnProtons) 
      _setupString.append("returnProtons 1 ");
    else
      _setupString.append("returnProtons 0 ");

    if (_returnGammas) 
      _setupString.append("returnGammas 1 ");
    else
      _setupString.append("returnGammas 0 ");

    if (_returnElectrons) 
      _setupString.append("returnElectrons 1 ");
    else
      _setupString.append("returnElectrons 0 ");

    if (_returnPions) 
      _setupString.append("returnPions 1 ");
    else
      _setupString.append("returnPions 0 ");

    if (_returnKaons) 
      _setupString.append("returnKaons 1 ");
    else
      _setupString.append("returnKaons 0 ");

    char tmpStr[256];
    sprintf(tmpStr, "date %d-%d-%d ", _month, _day, _year);
    _setupString.append(tmpStr);

    sprintf(tmpStr, "latitude %f ", _latitude);
    _setupString.append(tmpStr);
    
    sprintf(tmpStr, "altitude %d ", _altitude);
    _setupString.append(tmpStr);

    sprintf(tmpStr, "subboxLength %f ", _subboxLength);
    _setupString.append(tmpStr);

    _cryDataPath = std::string(std::getenv("CRYDATAPATH"));
    if( _cryDataPath.length() == 0) 
    {
      mf::LogError("CosmicCRY") << "no variable CRYDATAPATH set. Exit now.";
      exit(0);
    }

    _crySetup = new CRYSetup(_setupString, _cryDataPath);

    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(),
        &CLHEP::HepRandomEngine::flat);
    _crySetup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);

    if (_verbose > 1) 
      CLHEP::HepRandom::getTheEngine()->showStatus();

    _cryGen = new CRYGenerator(_crySetup);

    // Histograming
    if (_doHistograms) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("CosmicCRY");
      _hStartXZ = tfdir.make<TH2D>("StartXZ", "StartXZ", 
          500, -1.0e5,  1.0e5, 500, -1.0e5, 1.0e5 );
      _hStartY = tfdir.make<TH1D>("StartY", "StartY", 100, -1.0e3, 1.0e3 );

      // _hStartPlane = tfdir.make<TH1D>("StartPlane", "StartPlane", 5, 0, 5);
      _hStartE = tfdir.make<TH1D>("StartE", "StartE", 2000, 0, _muEMax);
      _hStartTheta = tfdir.make<TH1D>("StartTheta", "StartTheta", 200, -M_PI - 0.5, M_PI + 0.5);
      _hStartPhi = tfdir.make<TH1D>("StartPhi", "StartPhi", 200, -M_PI - 0.5, M_PI + 0.5);

    // And tree
      if (_saveTree) {
        _tCry = tfdir.make<TTree>("cryTree", "cryTree");
        _tCry->Branch("pdgId", &pdgId);
        _tCry->Branch("ke", &ke0);
        _tCry->Branch("px", &px0);
        _tCry->Branch("py", &py0);
        _tCry->Branch("pz", &pz0);
        _tCry->Branch("ptot", &ptot0);
        _tCry->Branch("x", &x0);
        _tCry->Branch("y", &y0);
        _tCry->Branch("z", &z0);
        _tCry->Branch("theta", &theta0);
        _tCry->Branch("phi", &phi0);
      }
    }
  }

  CosmicCRY::~CosmicCRY() { }

  void CosmicCRY::generate( GenParticleCollection& genParts )
  {
    std::vector<CRYParticle*> *daughters = new std::vector<CRYParticle*>;
    _cryGen->genEvent(daughters);
    for (unsigned j=0; j<daughters->size(); j++) {
      CRYParticle* p = (*daughters)[j];

      GlobalConstantsHandle<ParticleDataTable> pdt;
      const HepPDT::ParticleData& p_data = pdt->particle(p->PDGid()).ref();
      double mass = p_data.mass().value(); // in MeV

      double ke = p->ke(); // MeV by default in CRY
      double totalE = ke + mass;
      double totalP = safeSqrt(totalE * totalE - mass * mass);

      // Change coordinate system since y points upward, z points along
      // the beam line; which make cry(xyz) -> mu2e(zxy), uvw -> mu2e(zxy)
      CLHEP::Hep3Vector position(p->y() * 1000, p->z() * 1000, p->x() * 1000); // to mm
      CLHEP::HepLorentzVector mom4(totalP*p->v(), totalP*p->w(),
          totalP*p->u(), totalE);
      genParts.push_back(
          GenParticle(static_cast<PDGCode::type>(p->PDGid()), GenId::cosmicCRY,
            position, mom4, p->t()));

      if (_doHistograms) {
        _hStartXZ->Fill(position.x(), position.z());
        _hStartY->Fill(position.y());
        _hStartE->Fill(p->ke());
        // _hStartPlane->Fill(position.y());
        CLHEP::Hep3Vector momDir(p->u(), p->v(), p->w());
        _hStartTheta->Fill(momDir.theta());
        _hStartPhi->Fill(momDir.phi());

        pdgId = p->PDGid();
        ke0 = p->ke();
        x0 = p->x() * 1000;
        y0 = p->y() * 1000;
        z0 = p->z() * 1000;
        px0 = totalP * p->u();
        py0 = totalP * p->v();
        pz0 = totalP * p->w();
        ptot0 = totalP;
        theta0 = momDir.theta();
        phi0 = momDir.phi();

        _tCry->Fill();

      }

      if (_verbose > 1) {
        std::cout << "Secondary " << j 
          << ": " << CRYUtils::partName(p->id()) 
          << " (pdgId " << p->PDGid() << ")"
          << ", kinetic energy " << p->ke()  << " MeV"
          << ", position " 
          << "(" << p->x()
          << ", " << p->y()
          << ", " << p->z()
          << ") m"
          << ", time " << p->t() << " sec"
          << ", mass: " << mass
          << ", mom: " << totalP
          << ", mom dir.: " << p->u() <<", " << p->v() << ", " << p->w()
          << std::endl;

        // std::cout <<  genParts.back() << std::endl;
      }
    }
  }
}
