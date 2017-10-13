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
  , _refPointChoice(UNDEFINED)
  , _directionChoice(ALL)
  , _cosmicReferencePointInMu2e()
  , _vertical(false)
  , _dontProjectToSurface(config.getBool("cosmicCRY.dontProjectToSurface", false))

  {
    mf::LogInfo log("CosmicCRY");

    _cryDataPath = std::string(std::getenv("CRYDATAPATH"));
    if( _cryDataPath.length() == 0) 
    {
      mf::LogError("CosmicCRY") << "no variable CRYDATAPATH set. Exit now.";
      exit(0);
    }

    createSetupString();
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

    }

    // And tree
    if (_saveTree) {
      makeTrees();
    }
  }


  CosmicCRY::~CosmicCRY() { }

  void CosmicCRY::generate( GenParticleCollection& genParts )
  {
    std::vector<CRYParticle*> *daughters = new std::vector<CRYParticle*>;
    _cryGen->genEvent(daughters);

    if (_saveTree) {
      _pdgId0 = _cryGen->primaryParticle()->PDGid();
      _ke0 = _cryGen->primaryParticle()->ke(); // MeV
      _nSecondaries = daughters->size();
      _tCryPrimary->Fill();
    }

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
          GenParticle(static_cast<PDGCode::type>(p->PDGid()), GenId::cosmicDYB,
            position, mom4, p->t()));

      if (_doHistograms) {
        _hStartXZ->Fill(position.x(), position.z());
        _hStartY->Fill(position.y());
        _hStartE->Fill(p->ke());
        // _hStartPlane->Fill(position.y());
        // Keep the original uvw order to have correct theta and phi.
        CLHEP::Hep3Vector momDir(p->u(), p->v(), p->w());
        _hStartTheta->Fill(momDir.theta());
        _hStartPhi->Fill(momDir.phi());
      }

      if (_saveTree && j < _maxNSecondaries) {
        _pdgId1[j] = p->PDGid();
        _x1[j] = position.x();
        _y1[j] = position.y();
        _z1[j] = position.z();
        _t1[j] = p->t();
        _ke1[j] = p->ke();
        _px1[j] = totalP * p->v();
        _py1[j] = totalP * p->w();
        _pz1[j] = totalP * p->u();
        _ptot1[j] = totalP;

        // Keep the original uvw order to have correct theta and phi.
        CLHEP::Hep3Vector momDir(p->u(), p->v(), p->w()); 
        _theta1[j] = momDir.theta();
        _phi1[j] = momDir.phi();
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

    if (_saveTree) {
      _tCrySecondaries->Fill();
    }
  }

  void CosmicCRY::makeTrees()
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CosmicCRY");

    _tCryPrimary = tfdir.make<TTree>("cryPrimaryTree", "cryPrimaryTree");
    _tCryPrimary->Branch("evtId", &_evtId0);
    _tCryPrimary->Branch("pdgId", &_pdgId0);
    _tCryPrimary->Branch("ke", &_ke0);
    _tCryPrimary->Branch("nSecondaries", &_nSecondaries);

    _tCrySecondaries = tfdir.make<TTree>("crySecondariesTree",
        "crySecondariesTree");
    _tCrySecondaries->Branch("evtId", &_evtId0);
    _tCrySecondaries->Branch("nSecondaries", &_nSecondaries);
    _tCrySecondaries->Branch("pdgId", _pdgId1, "pdgId[nSecondaries]/I");
    _tCrySecondaries->Branch("x", _x1, "x[nSecondaries]/D");
    _tCrySecondaries->Branch("y", _y1, "y[nSecondaries]/D");
    _tCrySecondaries->Branch("z", _z1, "z[nSecondaries]/D");
    _tCrySecondaries->Branch("t", _t1, "t[nSecondaries]/D");
    _tCrySecondaries->Branch("ke", _ke1, "ke[nSecondaries]/D");
    _tCrySecondaries->Branch("px", _px1, "px[nSecondaries]/D");
    _tCrySecondaries->Branch("py", _py1, "py[nSecondaries]/D");
    _tCrySecondaries->Branch("pz", _pz1, "pz[nSecondaries]/D");
    _tCrySecondaries->Branch("ptot", _ptot1, "ptot[nSecondaries]/D");
    _tCrySecondaries->Branch("theta", _theta1, "theta[nSecondaries]/D");
    _tCrySecondaries->Branch("phi", _phi1, "phi[nSecondaries]/D");
  }

  void CosmicCRY::createSetupString()
  {
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
  }
}
